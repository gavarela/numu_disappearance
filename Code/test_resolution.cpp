#include <vector>

void test_resolution() {
    
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
    
    //// Take resolution into account
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Smear with a Gaussian and get covariance
    TF1 mygaus("mygaus", "[0] / sqrt(2*pi*[1]^2) * exp(-(x - [2])^2 / (2 * [1]^2))", 0, 10);
        // [0] is counts, [1] is counts*res and [2] energy
    
    // Get area - to normalise
    double area = 0;
    for (int i = 0; i < n; i++) {
        area += binwidth[i] * flux[i];
    }
    
    double resolution, newarea = 0;
    vector <double> resolutions = {0.05, 0.1, 0.25}, nzeros(n, 0);
    vector <vector <double> > fluxes(resolutions.size(), nzeros);
    for (int i = 0; i < resolutions.size(); i++) {
        resolution = resolutions[i];
        
        // Assign values (convolution)
        for (int k = 0; k < n; k++) {
            if (flux[k] != 0) {
                mygaus.SetParameters(flux[k], resolution*TMath::Sqrt(energy[k]), energy[k]);
                for (int j = 0; j < n; j++) {
                    fluxes[i][j] += mygaus(energy[j]);
                }
            }
        }
        
        // Get new area
        newarea = 0;
        for (int j = 0; j < n; j++) {
            newarea += fluxes[i][j] * binwidth[j];
        }
        
        // Normalise areas
        for (int j = 0; j < n; j++) {
            fluxes[i][j] = fluxes[i][j] / newarea * area;
        }
        
    }
    
    // Plot
    TCanvas *res_canvas = new TCanvas();
    
	TGraph *gr = new TGraph();
    for (int i = 0; i < n; i++) {
        gr->SetPoint(i, energy[i], flux[i]);
    }
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(30);
    
    TGraph *gr_05 = new TGraph();
    for (int i = 0; i < n; i++) {
        gr_05->SetPoint(i, energy[i], fluxes[0][i]);
    }
    gr_05->SetMarkerStyle(20);
    gr_05->SetMarkerSize(0.5);
    gr_05->SetMarkerColor(38);
    
    TGraph *gr_10 = new TGraph();
    for (int i = 0; i < n; i++) {
        gr_10->SetPoint(i, energy[i], fluxes[1][i]);
    }
    gr_10->SetMarkerStyle(20);
    gr_10->SetMarkerSize(0.5);
    gr_10->SetMarkerColor(46);
    
    TGraph *gr_25 = new TGraph();
    for (int i = 0; i < n; i++) {
        gr_25->SetPoint(i, energy[i], fluxes[2][i]);
    }
    gr_25->SetMarkerStyle(20);
    gr_25->SetMarkerSize(0.5);
    gr_25->SetMarkerColor(42);
    
    gr->SetTitle("Flux - Imperfect Resolutions; Energy (GeV); Flux (#Phi / m^{2} / 10^{6} POT)");
                    
    TLegend *res_legend = new TLegend();
    res_legend->AddEntry(gr, "Original", "p");
    res_legend->AddEntry(gr_05, "5%", "p");
    res_legend->AddEntry(gr_10, "10%", "p");
    res_legend->AddEntry(gr_25, "25%", "p");
    
    gr->Draw("AP");
    gr_05->Draw("P same");
    gr_10->Draw("P same");
    gr_25->Draw("P same");
    res_legend->Draw();
    
    res_canvas->SaveAs("../Out/Res/Flux_Res.root");
    res_canvas->SaveAs("../Out/Res/Flux_Res.png");
    
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

    double exposure = 6e20 /* POT */;
    
    // Initialise stuff
    double minchisq = 1e99, target = 23.40, np = 500;
    vector <double> flux0(n), xsec0(n), nearflux_null(n), neardetected_null(n), farflux_null(n), fardetected_null(n), dm2(np), sin2theta(np), corners(4, 0), twozeros(2, 0), npzeros(np, 0);
    vector <vector <double> > chisq(np, npzeros), chisq_diffs(np, npzeros), contour, sin_contour(resolutions.size()+1), dm2_contour(resolutions.size()+1), vecof2vecs(np, twozeros);
    vector <vector <vector <double> > > box_minmax(np, vecof2vecs);
    
    for (int i = 0; i < np; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(np-1));
        sin2theta[i] = TMath::Power(10, -3.0 + i*3.0/(np-1));
    }
    
    // Make prediction and get chisq and contours
    for  (int r = 0; r < resolutions.size()+1; r++) {
        
        // Null case - no oscillations
        for (int i = 0; i < n; i++) {
            xsec0[i] = xsec_forflux[i] * 1e-42;
            if (r != 3) {
                flux0[i] = fluxes[r][i] * (exposure / 1e6) * TMath::Power(MicroBooNE_dist, 2) / (1 - xsec0[i]*nAr_SBN);
            } else {
                flux0[i] = flux[i] * (exposure / 1e6) * TMath::Power(MicroBooNE_dist, 2) / (1 - xsec0[i]*nAr_SBN);
            }
            
            nearflux_null[i] = flux0[i] / TMath::Power(SBN_dist, 2);
            neardetected_null[i] = nearflux_null[i] * (xsec0[i]*nAr_SBN);
            
            farflux_null[i] = flux0[i] * (1 - xsec0[i]*nAr_SBN) / TMath::Power(ICARUS_dist, 2);
            fardetected_null[i] = farflux_null[i] * (xsec0[i]*nAr_ICARUS);
        }
        
        // Loop over phase space calculating Chisq
        clock_t startchi = clock();
        if (r != 3) {
            cout << "Getting chi squareds for resolution " << resolutions[r] << "." << endl;
        } else {
            cout << "Getting chi squareds for perfect resolution." << endl;
        }
        
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
        if (r != 3) {
            cout << "Getting contours for resolution " << resolutions[r] << "." << endl;
        } else {
            cout << "Getting contours for perfect resolution." << endl;
        }
        
        
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
    
    TGraph *gr_cont = new TGraph();
    for (int i = 0; i < sin_contour[3].size(); i++) {
        gr_cont->SetPoint(i, sin_contour[3][i], dm2_contour[3][i]);
    }
    gr_cont->SetMarkerStyle(20);
    gr_cont->SetMarkerSize(0.1);
    gr_cont->SetMarkerColor(30);
    gr_cont->SetLineColor(30);
    
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
    
    gr_cont->SetTitle("ICARUS 5#sigma CL - Imperfect Resolution; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend = new TLegend();
    legend->AddEntry(gr_cont, "No smearing", "l");
    legend->AddEntry(gr_05_cont, "5% smearing", "l");
    legend->AddEntry(gr_10_cont, "10% smearing", "l");
    legend->AddEntry(gr_25_cont, "25% smearing", "l");
    
    canvas->SetLogy();
    canvas->SetLogx();
    
    gr_cont->Draw("AP");
    gr_cont->GetXaxis()->SetLimits(0.001, 1);
    gr_cont->GetYaxis()->SetRangeUser(0.01, 100);
    canvas->Update();
    
    gr_05_cont->Draw("P same");
    gr_10_cont->Draw("P same");
    gr_25_cont->Draw("P same");
    legend->Draw();
    
    canvas->SaveAs("../Out/Res/Sensitivity_res.root");
    canvas->SaveAs("../Out/Res/Sensitivity_res.png");
    
    //// End timer
    //// ~~~~~~~~~
    
    clock_t end = clock();
    clock_t ticks = end - start;                        // in n of ticks
    double time = ticks / (double) CLOCKS_PER_SEC;      // make into secs
    
    cout << "Total time running: " << time << endl;
    
}