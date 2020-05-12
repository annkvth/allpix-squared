/**
 * @file
 * @brief Implementation of charge sensitive amplifier digitization module
 * @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "CSADigitizerModule.hpp"

#include "core/utils/unit.h"
#include "tools/ROOT.h"

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TProfile.h>

#include "objects/PixelHit.hpp"

using namespace allpix;

CSADigitizerModule::CSADigitizerModule(Configuration& config,
                                               Messenger* messenger,
                                               std::shared_ptr<Detector> detector)
    : Module(config, std::move(detector)), messenger_(messenger), pixel_message_(nullptr) {
    // Enable parallelization of this module if multithreading is enabled
    enable_parallelization();

    // Require PixelCharge message for single detector
    messenger_->bindSingle(this, &CSADigitizerModule::pixel_message_, MsgFlags::REQUIRED);

    // Seed the random generator with the global seed
    random_generator_.seed(getRandomSeed());

    // Set defaults for config variables
    config_.setDefault<double>("krummenacher_current", Units::get(20e-9, "C/s"));
    config_.setDefault<double>("feedback_capacitance", Units::get(5e-15, "C/V"));
    config_.setDefault<double>("feedback_resistance", Units::get(5e6, "V*s/C"));

    config_.setDefault<double>("parasitic_capacitance", 50.0);
    config_.setDefault<double>("transconductance", 1.0);

    config_.setDefault<int>("threshold", Units::get(800, "mV"));
    config_.setDefault<int>("threshold_smearing", Units::get(5, "mV"));

    config_.setDefault<int>("clock_cycle", Units::get(10, "ns"));

    config_.setDefault<bool>("output_pulsegraphs", false);
    config_.setDefault<bool>("output_plots", config_.get<bool>("output_pulsegraphs"));
    config_.setDefault<int>("output_plots_scale", Units::get(30, "ke"));
    config_.setDefault<int>("output_plots_bins", 100);


    // Copy some variables from configuration to avoid lookups:
    cf_ = config_.get<double>("feedback_capacitance");
    ct_ = config_.get<double>("parasitic_capacitance");
    g_ = config_.get<double>("transconductance");
    ikrum_ = config_.get<double>("krummenacher_current");

    
    output_plots_ = config_.get<bool>("output_plots");
    output_pulsegraphs_ = config_.get<bool>("output_pulsegraphs");

}


//asv init just a copy of generic digitizer for now - histograms don't make too much sense
void CSADigitizerModule::init() {

    if(config_.get<bool>("output_plots")) {
        LOG(TRACE) << "Creating output plots";

        // Plot axis are in kilo electrons - convert from framework units!
        int maximum = static_cast<int>(Units::convert(config_.get<int>("output_plots_scale"), "ke"));
        auto nbins = config_.get<int>("output_plots_bins");

        // Create histograms if needed
        h_pxq = new TH1D("pixelcharge", "raw pixel charge;pixel charge [ke];pixels", nbins, 0, maximum);
	h_amplified_charge_ = new TH1D("amplifiedcharge", "amplified charge;amplified pixel charge [ke];events", nbins, 0, maximum);
        h_pulse_charge_ = new TH1D("pulsecharge", "input pulse charge per pixel;input pulse pixel charge [ke];pixels", nbins, 0, maximum);

    }
}

void CSADigitizerModule::run(unsigned int event_num) {
    // Loop through all pixels with charges
    std::vector<PixelHit> hits;
    for(auto& pixel_charge : pixel_message_->getData()) {
        auto pixel = pixel_charge.getPixel();
        auto pixel_index = pixel.getIndex();
        auto inputcharge = static_cast<double>(pixel_charge.getCharge());

        LOG(DEBUG) << "Received pixel " << pixel_index << ", charge " << Units::display(inputcharge, "e");

	//asv copied from Simon:
	const auto& pulse = pixel_charge.getPulse();
	auto pulse_vec = pulse.getPulse();
	
        auto transfer = [&](double q_ind, double time) {
            // Transfer function for Krummenacher circuit:
            // h(t) = Q / C_f * exp (w_2 * t − exp(w_1*t))
            // with w_1 = g * C_f / C_t and w_2 = 1 / C_f / R_f = 1 / C_f / I_krum * 20
            return (q_ind / cf_ * exp(time / cf_ / ikrum_ * 20 - exp(g_ * cf_ / ct_ * time)));
        };
	

        double charge = 0;
        size_t steps = 0;
        auto timestep = pulse.getBinning();
	LOG(TRACE) << "Preparing pulse for pixel " << pixel_index << ", " << pulse_vec.size() << " bins of "
                       << Units::display(timestep, {"ps", "ns"})
                       << ", total charge: " << Units::display(pulse.getCharge(), "e");

	
        for(auto& q_ind : pulse_vec) {
            LOG(TRACE) << "Charge " << Units::display(charge, "e") << " Transfer at "
                       << timestep * static_cast<double>(steps) << ": "
                       << transfer(charge, timestep * static_cast<double>(steps)) << " "
                       << transfer(q_ind, timestep * static_cast<double>(steps));
            charge += q_ind;
            charge -= ikrum_ * timestep;
            steps++;
        }
	
        if(config_.get<bool>("output_plots")) {
            h_pxq->Fill(inputcharge / 1e3);
	    h_pulse_charge_->Fill(pulse.getCharge() / 1e3);
	    h_amplified_charge_->Fill(charge / 1e3);
        }

        // Fill a graphs with the individual pixel pulses:
        if(output_pulsegraphs_) {
            // Generate x-axis:
            std::vector<double> time(pulse_vec.size());
            // clang-format off
            std::generate(time.begin(), time.end(), [n = 0.0, timestep]() mutable {  auto now = n; n += timestep; return now; });
            // clang-format on

            std::string name =
                "pulse_ev" + std::to_string(event_num) + "_px" + std::to_string(pixel_index.x()) + "-" + std::to_string(pixel_index.y());
            auto pulse_graph = new TGraph(static_cast<int>(pulse_vec.size()), &time[0], &pulse_vec[0]);
            pulse_graph->GetXaxis()->SetTitle("t [ns]");
            pulse_graph->GetYaxis()->SetTitle("Q_{ind} [e]");
            pulse_graph->SetTitle(("Induced charge in pixel (" + std::to_string(pixel_index.x()) + "," +
                                   std::to_string(pixel_index.y()) +
                                   "), Q_{tot} = " + std::to_string(pulse.getCharge()) + " e")
                                      .c_str());
            getROOTDirectory()->WriteTObject(pulse_graph, name.c_str());

            // Generate graphs of integrated charge over time:
            std::vector<double> charge_vec;
            charge = 0;
            for(const auto& bin : pulse_vec) {
                charge += bin;
                charge_vec.push_back(charge);
            }

            name = "charge_ev" + std::to_string(event_num) + "_px" + std::to_string(pixel_index.x()) + "-" +
                   std::to_string(pixel_index.y());
            auto charge_graph = new TGraph(static_cast<int>(charge_vec.size()), &time[0], &charge_vec[0]);
            charge_graph->GetXaxis()->SetTitle("t [ns]");
            charge_graph->GetYaxis()->SetTitle("Q_{tot} [e]");
            charge_graph->SetTitle(("Accumulated induced charge in pixel (" + std::to_string(pixel_index.x()) + "," +
                                    std::to_string(pixel_index.y()) +
                                    "), Q_{tot} = " + std::to_string(pulse.getCharge()) + " e")
                                       .c_str());
            getROOTDirectory()->WriteTObject(charge_graph, name.c_str());
        }


	
        // Add the hit to the hitmap
        hits.emplace_back(pixel, 0, charge, &pixel_charge);
    }

    // Output summary and update statistics
    LOG(INFO) << "Digitized " << hits.size() << " pixel hits";
    total_hits_ += hits.size();

    if(!hits.empty()) {
        // Create and dispatch hit message
        auto hits_message = std::make_shared<PixelHitMessage>(std::move(hits), getDetector());
        messenger_->dispatchMessage(this, hits_message);
    }
}

void CSADigitizerModule::finalize() {
    if(config_.get<bool>("output_plots")) {
        // Write histograms
        LOG(TRACE) << "Writing output plots to file";
        h_pxq->Write();
        h_pulse_charge_->Write();
        h_amplified_charge_->Write();
	
        // if(config_.get<int>("adc_resolution") > 0) {
        //     h_pxq_adc_smear->Write();
        //     h_calibration->Write();
        // }
    }

    LOG(INFO) << "Digitized " << total_hits_ << " pixel hits in total";
}
