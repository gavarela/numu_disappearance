#include <vector>

void test_POTs() {
    
    //// Start Timer
    //// ~~~~~~~~~~~
    
    clock_t start = clock(); // n of ticks since start of program
    
    //// Get Data
    //// ~~~~~~~~
    
    // From flux histogram
    TFile fluxfile("../In/flux.root");
    TH1D *fluxhist = (TH1D*)fluxfile.Get("numu_CV_AV_TPC");

    int n = fluxhist->GetSize() - 2;
    vector <double> energy(n), flux(n), binmin(n), binwidth(n);
    for (int i = 0; i < n; i++) {
        energy[i] = fluxhist->GetBinCenter(i+1);
        flux[i] = fluxhist->GetBinContent(i+1);
        binmin[i] = fluxhist->GetBinLowEdge(i+1); // All 0.025 less than energy
        binwidth[i] = fluxhist->GetBinWidth(i+1); // All 0.05
    }
    
    // From cross-section graph
    TFile xsecfile("../In/cross_section.root");
    TGraph *xsecgraph = (TGraph*)xsecfile.Get("qel_cc_n");
    
    int m = xsecgraph->GetN();
    vector <double> energyx(m), xsec(m);
    for (int i = 0; i < m; i++) {
        xsecgraph->GetPoint(i, energyx[i], xsec[i]);
    }
    
    // Interpolate to make xsec the same length as flux
    vector <double> xsec_forflux(n, 0);
    int energyx_ind = 0;
    for (int i = 0; i < n; i++) {
        while (energy[i] > energyx[energyx_ind]) {
            energyx_ind += 1;
        }
        xsec_forflux[i] = xsec[energyx_ind-1] + (xsec[energyx_ind] - xsec[energyx_ind-1])*(energy[i] - energyx[energyx_ind-1])/(energyx[energyx_ind] - energyx[energyx_ind-1]);
    }
    
    //// Get 5sigma contour - skip plots this time
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Oscillation function
    TF1 nutonu("nutonu", "1 - [0] * (sin(1.27 * [1] * x))^2", 0, 25); // [0] is sin^2(2theta), [1] is delta m^2, x is L/E
    
    // Some parameters:
    double SBN_dist = 0.1, MicroBooNE_dist = 0.47, ICARUS_dist = 0.6; // all in km
    
    double ICARUS_volume = 3.2 * 2.96 * 18 * 2 /* m^3 */, 
        LAr_density = 1400 /* kg / m^3 */,
        Ar_atomicmass = 39.948 * 1.66054e-27 /* kg */,
        nAr_ICARUS = ICARUS_volume * LAr_density / Ar_atomicmass,
        SBN_volume = 4 * 4 * 5 /* m^3 */,
        nAr_SBN = SBN_volume * LAr_density / Ar_atomicmass;
    
    // Initialise stuff
    vector <double> exposures = {6e22, 6e20, 6e18};
    
    double minchisq = 1e99, target = 23.40, np = 500, exposure;
    vector <double> flux0(n), xsec0(n), nearflux_null(n), neardetected_null(n), farflux_null(n), fardetected_null(n), dm2(np), sin2theta(np), corners(4, 0), twozeros(2, 0), npzeros(np, 0);
    vector <vector <double> > chisq(np, npzeros), chisq_diffs(np, npzeros), contour, sin_contour(exposures.size()), dm2_contour(exposures.size()), vecof2vecs(np, twozeros);
    vector <vector <vector <double> > > box_minmax(np, vecof2vecs);
    
    for (int i = 0; i < np; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(np-1));
        sin2theta[i] = TMath::Power(10, -4.0 + i*4.0/(np-1));
    }
    
    // Make prediction and get chisq and contours
    for  (int r = 0; r < exposures.size(); r++) {
        
        exposure = exposures[r];
        
        // Null case - no oscillations
        for (int i = 0; i < n; i++) {
            xsec0[i] = xsec_forflux[i] * 1e-42;
            flux0[i] = flux[i] * (exposure / 1e6) * TMath::Power(MicroBooNE_dist, 2) / (1 - xsec0[i]*nAr_SBN);
        
            nearflux_null[i] = flux0[i] / TMath::Power(SBN_dist, 2);
            neardetected_null[i] = nearflux_null[i] * (xsec0[i]*nAr_SBN);
            
            farflux_null[i] = flux0[i] * (1 - xsec0[i]*nAr_SBN) / TMath::Power(ICARUS_dist, 2);
            fardetected_null[i] = farflux_null[i] * (xsec0[i]*nAr_ICARUS);
        }
        
        // Loop over phase space calculating Chisq
        clock_t startchi = clock();
        cout << "Getting chi squareds for exposure " << exposure << "." << endl;
        
        minchisq = 1e99;
        for (int i = 0; i < np; i++){
            for (int j = 0; j < np; j++) {
                
                // Set function parameters
                nutonu.SetParameters(sin2theta[i], dm2[j]);
                
                // Find oscillation flux and detections and calculate chisq
                chisq[i][j] = 0;
                
                for (int k = 0; k < n; k++) {
                    if (neardetected_null[k] != 0) {
                        chisq[i][j] += fardetected_null[k] * TMath::Power(1 - nutonu(ICARUS_dist/energy[k]), 2) / (1 + fardetected_null[k]/neardetected_null[k]);
                    }
                }
                
                if (chisq[i][j] < minchisq) {
                    minchisq = chisq[i][j];
                }
            
            }
        }
    
        clock_t endchi = clock();
        clock_t tickschi = endchi - startchi;                    // in n of ticks
        double timechi = tickschi / (double) CLOCKS_PER_SEC;     // make into secs
        
        cout << "Done in " << timechi << "s. " << endl;
        
        // Get differences
        for (int i = 0; i < chisq.size(); i++) {
            for (int j = 0; j < chisq[0].size(); j++) {
                chisq_diffs[i][j] = chisq[i][j] - minchisq;
            }
        }
        
        // Get contours
        clock_t startcont = clock();
        cout << "Getting contours for exposure " << exposure << "." << endl;
        
        // Initialise box_minmax
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                box_minmax[i][j] = {0, 0};
            }
        }
        
        // Fill box_minmax out
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                
                corners = {chisq_diffs[i][j], chisq_diffs[i+1][j], chisq_diffs[i][j+1], chisq_diffs[i+1][j+1]};
                box_minmax[i][j] = {corners[0], corners[0]};
                for (int l = 1; l < 4; l++) {
                    if (corners[l] > box_minmax[i][j][1]) {
                        box_minmax[i][j][1] = corners[l];
                    } else if (corners[l] < box_minmax[i][j][0]) {
                        box_minmax[i][j][0] = corners[l];
                    }
                }
                
            }
        }
        
        // Get the contour
        contour.clear();
        for (int i = 0; i < np-1; i++) {
            for (int j = 0; j < np-1; j++) {
                
                if ((target >= box_minmax[i][j][0]) && (target <= box_minmax[i][j][1])) {
                    contour.push_back({(sin2theta[i] + sin2theta[i+1])/2, (dm2[j] + dm2[j+1])/2});
                }
            
            }
        }
        
        // Save the contour
        for (int j = 0; j < contour.size(); j++) {
            sin_contour[r].push_back(contour[j][0]);
            dm2_contour[r].push_back(contour[j][1]);
        }
        
        clock_t endcont = clock();
        clock_t tickscont = endcont - startcont;                    // in n of ticks
        double timecont = tickscont / (double) CLOCKS_PER_SEC;      // make into secs
        
        cout << "Done in " << timecont << "s." << endl;
        
        cout << "Sin contour sizes: " << sin_contour[0].size() << ", " << sin_contour[1].size() << ", " << sin_contour[2].size() << endl;
        cout << "dm2 contour sizes: " << dm2_contour[0].size() << ", " << dm2_contour[1].size() << ", " << dm2_contour[2].size() << endl;
        
    }
    
    // Now plot it
    TCanvas *canvas = new TCanvas();
    
	TGraph *gr_05_cont = new TGraph();
    for (int i = 0; i < sin_contour[0].size(); i++) {
        gr_05_cont->SetPoint(i, sin_contour[0][i], dm2_contour[0][i]);
    }
    gr_05_cont->SetMarkerStyle(20);
    gr_05_cont->SetMarkerSize(0.1);
    gr_05_cont->SetMarkerColor(38);
    gr_05_cont->SetLineColor(38);
    
    TGraph *gr_10_cont = new TGraph();
    for (int i = 0; i < sin_contour[1].size(); i++) {
        gr_10_cont->SetPoint(i, sin_contour[1][i], dm2_contour[1][i]);
    }
    gr_10_cont->SetMarkerStyle(20);
    gr_10_cont->SetMarkerSize(0.1);
    gr_10_cont->SetMarkerColor(46);
    gr_10_cont->SetLineColor(46);
    
    TGraph *gr_25_cont = new TGraph();
    for (int i = 0; i < sin_contour[2].size(); i++) {
        gr_25_cont->SetPoint(i, sin_contour[2][i], dm2_contour[2][i]);
    }
    gr_25_cont->SetMarkerStyle(20);
    gr_25_cont->SetMarkerSize(0.1);
    gr_25_cont->SetMarkerColor(42);
    gr_25_cont->SetLineColor(42);
    
    gr_05_cont->SetTitle("ICARUS 5#sigma CL - Varying Exposure; sin^{2}(2#theta); #Delta m^{2}");
    
    TLegend *legend = new TLegend();
    legend->AddEntry(gr_05_cont, "6x10^{22} POT", "l");
    legend->AddEntry(gr_10_cont, "6x10^{20} POT", "l");
    legend->AddEntry(gr_25_cont, "6x10^{18} POT", "l");
    
    canvas->SetLogy();
    canvas->SetLogx();
    
    gr_05_cont->Draw("AP");
    gr_05_cont->GetXaxis()->SetLimits(0.001, 1);
    gr_05_cont->GetYaxis()->SetRangeUser(0.01, 100);
    canvas->Update();
    
    gr_10_cont->Draw("P same");
    gr_25_cont->Draw("P same");
    legend->Draw();
    
    canvas->SaveAs("../Out/POTs/Sensitivity_POTs.root");
    canvas->SaveAs("../Out/POTs/Sensitivity_POTs.png");
    
    //// When does ICARUS > MicroBooNE?
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Parameters and phase space
    double MicroBooNE_volume = 2.33 * 2.56 * 10.37 /* m^3 */,
        nAr_MicroBooNE = MicroBooNE_volume * LAr_density / Ar_atomicmass;
    
    cout << nAr_MicroBooNE << "   " << nAr_ICARUS << endl;
    
    exposures.clear();
    int n_exp = 300;
    for (int i = 0; i < n_exp; i++) {
        exposures.push_back(TMath::Power(10, 18.0 + i*(25.0 - 10.0)/(n_exp-1)));
    }
    
    // Get null counts at MicroBooNE and ICARUS
    double dmedflux_null, dmeddetected_null, dfarflux_null, dfardetected_null;
    vector <double> meddetected_null;
    meddetected_null.clear(); fardetected_null.clear();
    for (int r = 0; r < exposures.size(); r++) {
        
        exposure = exposures[r];
        dmeddetected_null = 0; dfardetected_null = 0;
        
        // Null case - no oscillations
        for (int i = 0; i < n; i++) {
            xsec0[i] = xsec_forflux[i] * 1e-42;
            
            dmedflux_null = flux[i] * ((13.2e20 + exposure) / 1e6);
            dmeddetected_null += dmedflux_null * (xsec0[i]*nAr_MicroBooNE);
            
            dfarflux_null = flux[i] * (exposure / 1e6) * (1 - xsec0[i]*nAr_MicroBooNE) * TMath::Power(MicroBooNE_dist/ICARUS_dist, 2);
            dfardetected_null += dfarflux_null * (xsec0[i]*nAr_ICARUS);
        }
        
        meddetected_null.push_back(dmeddetected_null);
        fardetected_null.push_back(dfardetected_null);
        
    }
    
    // Plot
    TCanvas *canvas_2ds = new TCanvas();
    
	TGraph *gr_MBNE = new TGraph();
    for (int i = 0; i < exposures.size(); i++) {
        gr_MBNE->SetPoint(i, exposures[i], meddetected_null[i]);
    }
    gr_MBNE->SetMarkerStyle(20);
    gr_MBNE->SetMarkerSize(0.1);
    gr_MBNE->SetMarkerColor(38);
    gr_MBNE->SetLineColor(38);
    
    TGraph *gr_ICARUS = new TGraph();
    for (int i = 0; i < exposures.size(); i++) {
        gr_ICARUS->SetPoint(i, exposures[i], fardetected_null[i]);
    }
    gr_ICARUS->SetMarkerStyle(20);
    gr_ICARUS->SetMarkerSize(0.1);
    gr_ICARUS->SetMarkerColor(46);
    gr_ICARUS->SetLineColor(46);
    
    gr_ICARUS->SetTitle("Volume of Data - MicroBooNE and ICARUS; Exposure (POT); Counts, N_{i}^{null}");
    
    TLegend *legend_2ds = new TLegend();
    legend_2ds->AddEntry(gr_05_cont, "MicroBooNE", "l");
    legend_2ds->AddEntry(gr_10_cont, "ICARUS", "l");
    
    gr_ICARUS->Draw("AP");
    gr_MBNE->Draw("P same");
    
    canvas_2ds->SetLogx();
    canvas_2ds->SetLogy();
    legend_2ds->Draw();
    
    canvas_2ds->SaveAs("../Out/POTs/MBNE_ICARUS.root");
    canvas_2ds->SaveAs("../Out/POTs/MBNE_ICARUS.png");
    
    
    //// When does this cover the best fit point?
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    cout << "Finding exposures to cover best fit point." << endl;
    double dm2_bf = 1.7, sin_bf = 0.062;
    
    // Make prediction and get chisq and contours
    vector <int> sigmas = {3, 5};
    vector <vector <vector <double> > > sin_contours(2), dm2_contours(2);
    for (int s = 0; s < sigmas.size(); s++) {
        
        if (s == 0) { // 3 sigma
            target = 7.75;
            exposures = {6e20, 6e19, 6e18};
            sin_contours.push_back({}); dm2_contour.push_back({});
        } else if (s == 1) {
            target = 23.4;
            exposures = {6e20, 6e19, 1.7e19};
            sin_contours.push_back({}); dm2_contour.push_back({});
        }
        
        for (int r = 0; r < exposures.size(); r++) {
            
            exposure = exposures[r];

            // Null case - no oscillations
            for (int i = 0; i < n; i++) {
                xsec0[i] = xsec_forflux[i] * 1e-42;
                flux0[i] = flux[i] * (exposure / 1e6) * TMath::Power(MicroBooNE_dist, 2) / (1 - xsec0[i]*nAr_SBN);

                nearflux_null[i] = flux0[i] / TMath::Power(SBN_dist, 2);
                neardetected_null[i] = nearflux_null[i] * (xsec0[i]*nAr_SBN);

                farflux_null[i] = flux0[i] * (1 - xsec0[i]*nAr_SBN) / TMath::Power(ICARUS_dist, 2);
                fardetected_null[i] = farflux_null[i] * (xsec0[i]*nAr_ICARUS);
            }

            // Loop over phase space calculating Chisq
            clock_t startchi = clock();
            cout << "Getting chi squareds for CL " << sigmas[s] << " sigma and exposure " << exposure << "." << endl;

            minchisq = 1e99;
            for (int i = 0; i < np; i++){
                for (int j = 0; j < np; j++) {

                    // Set function parameters
                    nutonu.SetParameters(sin2theta[i], dm2[j]);

                    // Find oscillation flux and detections and calculate chisq
                    chisq[i][j] = 0;

                    for (int k = 0; k < n; k++) {
                        if (neardetected_null[k] != 0) {
                            chisq[i][j] += fardetected_null[k] * TMath::Power(1 - nutonu(ICARUS_dist/energy[k]), 2) / (1 + fardetected_null[k]/neardetected_null[k]);
                        }
                    }

                    if (chisq[i][j] < minchisq) {
                        minchisq = chisq[i][j];
                    }

                }
            }

            clock_t endchi = clock();
            clock_t tickschi = endchi - startchi;                    // in n of ticks
            double timechi = tickschi / (double) CLOCKS_PER_SEC;     // make into secs

            cout << "Done in " << timechi << "s. " << endl;

            // Get differences
            for (int i = 0; i < chisq.size(); i++) {
                for (int j = 0; j < chisq[0].size(); j++) {
                    chisq_diffs[i][j] = chisq[i][j] - minchisq;
                }
            }

            // Get contours
            clock_t startcont = clock();
            cout << "Getting contours for CL " << sigmas[s] << " sigma and exposure " << exposure << "." << endl;

            // Initialise box_minmax
            for (int i = 0; i < np-1; i++) {
                for (int j = 0; j < np-1; j++) {
                    box_minmax[i][j] = {0, 0};
                }
            }
            
            // Fill box_minmax out
            for (int i = 0; i < np-1; i++) {
                for (int j = 0; j < np-1; j++) {

                    corners = {chisq_diffs[i][j], chisq_diffs[i+1][j], chisq_diffs[i][j+1], chisq_diffs[i+1][j+1]};
                    box_minmax[i][j] = {corners[0], corners[0]};
                    for (int l = 1; l < 4; l++) {
                        if (corners[l] > box_minmax[i][j][1]) {
                            box_minmax[i][j][1] = corners[l];
                        } else if (corners[l] < box_minmax[i][j][0]) {
                            box_minmax[i][j][0] = corners[l];
                        }
                    }

                }
            }
            
            // Get the contour
            contour.clear();
            for (int i = 0; i < np-1; i++) {
                for (int j = 0; j < np-1; j++) {

                    if ((target >= box_minmax[i][j][0]) && (target <= box_minmax[i][j][1])) {
                        contour.push_back({(sin2theta[i] + sin2theta[i+1])/2, (dm2[j] + dm2[j+1])/2});
                    }

                }
            }
            
            // Save the contour
            sin_contours[s].push_back({}); dm2_contours[s].push_back({});
            for (int j = 0; j < contour.size(); j++) {
                sin_contours[s][r].push_back(contour[j][0]);
                dm2_contours[s][r].push_back(contour[j][1]);
            }
            
            clock_t endcont = clock();
            clock_t tickscont = endcont - startcont;                    // in n of ticks
            double timecont = tickscont / (double) CLOCKS_PER_SEC;      // make into secs

            cout << "Done in " << timecont << "s." << endl;
            
            
            cout << "For CL " << sigmas[s] << " sigma:" << endl;
            cout << "Sin contour sizes: " << sin_contours[s][0].size() << ", " << sin_contours[s][1].size() << ", " << sin_contours[s][2].size() << endl;
            cout << "dm2 contour sizes: " << dm2_contours[s][0].size() << ", " << dm2_contours[s][1].size() << ", " << dm2_contours[s][2].size() << endl;

        }
        
    }
    
    // Plot 3 sigma
    TCanvas *canvas_3s = new TCanvas();
    
	TGraph *gr_3s_0 = new TGraph();
    for (int i = 0; i < sin_contours[0][0].size(); i++) {
        gr_3s_0->SetPoint(i, sin_contours[0][0][i], dm2_contours[0][0][i]);
    }
    gr_3s_0->SetMarkerStyle(20);
    gr_3s_0->SetMarkerSize(0.1);
    gr_3s_0->SetMarkerColor(30);
    gr_3s_0->SetLineColor(30);
    
    TGraph *gr_3s_1 = new TGraph();
    for (int i = 0; i < sin_contours[0][1].size(); i++) {
        gr_3s_1->SetPoint(i, sin_contours[0][1][i], dm2_contours[0][1][i]);
    }
    gr_3s_1->SetMarkerStyle(20);
    gr_3s_1->SetMarkerSize(0.1);
    gr_3s_1->SetMarkerColor(38);
    gr_3s_1->SetLineColor(38);
    
    TGraph *gr_3s_2 = new TGraph();
    for (int i = 0; i < sin_contours[0][2].size(); i++) {
        gr_3s_2->SetPoint(i, sin_contours[0][2][i], dm2_contours[0][2][i]);
    }
    gr_3s_2->SetMarkerStyle(20);
    gr_3s_2->SetMarkerSize(0.1);
    gr_3s_2->SetMarkerColor(46);
    gr_3s_2->SetLineColor(46);
    
    TGraph *gr_bestfit = new TGraph();
    gr_bestfit->SetPoint(0, sin_bf, dm2_bf);
    gr_bestfit->SetMarkerStyle(3);
    gr_bestfit->SetMarkerSize(0.7);
    gr_bestfit->SetMarkerColor(1);
    
    gr_3s_0->SetTitle("Best Fit Point - ICARUS 3#sigma Contour; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend_3s = new TLegend();
    legend_3s->AddEntry(gr_3s_0, "6x10^{20} POT", "l");
    legend_3s->AddEntry(gr_3s_1, "6x10^{19} POT", "l");
    legend_3s->AddEntry(gr_3s_2, "6x10^{18} POT", "l");
    legend_3s->AddEntry(gr_bestfit, "Best Fit Point", "p");
    
    gr_3s_0->Draw("AP");
    gr_3s_1->Draw("P same");
    gr_3s_2->Draw("P same");
    gr_bestfit->Draw("P same");
    
    canvas_3s->SetLogx();
    canvas_3s->SetLogy();
    legend_3s->Draw();
    
    canvas_3s->SaveAs("../Out/POTs/3sigma_bf.root");
    canvas_3s->SaveAs("../Out/POTs/3sigma_bf.png");
    
    // Plot 5 sigma
    TCanvas *canvas_5s = new TCanvas();
    
	TGraph *gr_5s_0 = new TGraph();
    for (int i = 0; i < sin_contours[1][0].size(); i++) {
        gr_5s_0->SetPoint(i, sin_contours[1][0][i], dm2_contours[1][0][i]);
    }
    gr_5s_0->SetMarkerStyle(20);
    gr_5s_0->SetMarkerSize(0.1);
    gr_5s_0->SetMarkerColor(30);
    gr_5s_0->SetLineColor(30);
    
    TGraph *gr_5s_1 = new TGraph();
    for (int i = 0; i < sin_contours[1][1].size(); i++) {
        gr_5s_1->SetPoint(i, sin_contours[1][1][i], dm2_contours[1][1][i]);
    }
    gr_5s_1->SetMarkerStyle(20);
    gr_5s_1->SetMarkerSize(0.1);
    gr_5s_1->SetMarkerColor(38);
    gr_5s_1->SetLineColor(38);
    
    TGraph *gr_5s_2 = new TGraph();
    for (int i = 0; i < sin_contours[1][2].size(); i++) {
        gr_5s_2->SetPoint(i, sin_contours[1][2][i], dm2_contours[1][2][i]);
    }
    gr_5s_2->SetMarkerStyle(20);
    gr_5s_2->SetMarkerSize(0.1);
    gr_5s_2->SetMarkerColor(46);
    gr_5s_2->SetLineColor(46);
    
    gr_5s_0->SetTitle("Best Fit Point - ICARUS 5#sigma Contour; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend_5s = new TLegend();
    legend_5s->AddEntry(gr_5s_0, "6x10^{20} POT", "l");
    legend_5s->AddEntry(gr_5s_1, "6x10^{19} POT", "l");
    legend_5s->AddEntry(gr_5s_2, "1.7x10^{19} POT", "l");
    legend_5s->AddEntry(gr_bestfit, "Best Fit Point", "p");
    
    gr_5s_0->Draw("AP");
    gr_5s_1->Draw("P same");
    gr_5s_2->Draw("P same");
    gr_bestfit->Draw("P same");
    
    canvas_5s->SetLogx();
    canvas_5s->SetLogy();
    legend_5s->Draw();
    
    canvas_5s->SaveAs("../Out/POTs/5sigma_bf.root");
    canvas_5s->SaveAs("../Out/POTs/5sigma_bf.png");
    
    
    //// When does the ND constraint match 10%, etc.?
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Still need to do...
    
    
    //// End timer
    //// ~~~~~~~~~
    
    clock_t end = clock();
    clock_t ticks = end - start;                        // in n of ticks
    double time = ticks / (double) CLOCKS_PER_SEC;      // make into secs
    
    cout << "Total time running: " << time << endl;
    
}