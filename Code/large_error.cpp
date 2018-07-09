#include <vector>

void large_error() {
    
    /// Start Timer
    //// ~~~~~~~~~~
    
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
    
    //// Get chi squareds and contours
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    // Oscillation function
    TF1 nutonu("nutonu", "1 - [0] * (sin(1.27 * [1] * x))^2", 0, 25); // [0] is sin^2(2theta), [1] is delta m^2, x is L/E
    
    // Some parameters:
    double MicroBooNE_dist = 0.47, ICARUS_dist = 0.6; // all in km
    
    double ICARUS_volume = 3.2 * 2.96 * 18 * 2 /* m^3 */, 
        LAr_density = 1400 /* kg / m^3 */,
        Ar_atomicmass = 39.948 * 1.66054e-27 /* kg */,
        nAr_atoms = ICARUS_volume * LAr_density / Ar_atomicmass;

    double exposure = 6e20 /* POT */;
    
    // Null case - no oscillations
    vector <double> farflux_null(n), fardetected_null(n);
    for (int k = 0; k < n; k++) {
        farflux_null[k] = flux[k] * (exposure / 1e6) * TMath::Power(MicroBooNE_dist/ICARUS_dist, 2);
        fardetected_null[k] = farflux_null[k] * (xsec_forflux[k] * 1e-42) * nAr_atoms;
    }
    
    // Phase space
    int np = 500;
    vector <double> dm2(np), sin2theta(np);
    for (int i = 0; i < np; i++) {
        dm2[i] = TMath::Power(10, -2.0 + i*4.0/(np-1));
        sin2theta[i] = TMath::Power(10, -2.0 + i*2.0/(np-1));
    }
    
    // Loop over phase space calculating Chisq
    clock_t startchi = clock();
    cout << "Calculating chi squareds..." << endl;
    
    double minchisq = 1e99, fardetected_osc;
    vector <double> npzeros(np, 0);
    vector <vector <double> > chisq(np, npzeros);
    for (int i = 0; i < np; i++){
        for (int j = 0; j < np; j++) {
        
            // Set function parameters
            nutonu.SetParameters(sin2theta[i], dm2[j]);
            
            // Find null and oscillation fluxes and detections and calculate chisq
            for (int k = 0; k < n; k++) {
                chisq[i][j] += fardetected_null[k] * TMath::Power(1 - nutonu(ICARUS_dist/energy[k]), 2) / (1 + 0.09*fardetected_null[k]);
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
    
    // Plot
    TCanvas *chisqcanvas = new TCanvas();
    
    TGraph2D *logchisqplot = new TGraph2D();
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            logchisqplot->SetPoint(i*np + j, TMath::Log10(sin2theta[i]), TMath::Log10(dm2[j]), chisq[i][j]);
        }
    }
    
    logchisqplot->SetTitle("#chi^{2} - 30% Knowledge of Flux ; log_{10}(sin^{2}(2#theta)); log_{10}(#Delta m^{2}); #chi^{2}");
    gStyle->SetPalette(1);
    logchisqplot->Draw("surf1");
    chisqcanvas->SaveAs("../Out/LargeError/chisqplot_largeerrors.root");
    
    // Get differences
    vector <vector <double> > chisq_diffs(np, npzeros);
    for (int i = 0; i < chisq.size(); i++) {
        for (int j = 0; j < chisq[0].size(); j++) {
            chisq_diffs[i][j] = chisq[i][j] - minchisq;
        }
    }
    
    // Get contours
    clock_t startcont = clock();
    cout << "Getting contours." << endl;
    
    double target;
    vector <double> target_dchisq = {1.64, 7.75, 23.40}, corners(4, 0), twozeros(2, 0);
    vector <vector <double> > contour, sin_contour(target_dchisq.size()), dm2_contour(target_dchisq.size()), vecof2vecs(np-1, twozeros);
    vector <vector <vector <double> > > box_minmax(np-1, vecof2vecs);
    for  (int k = 0; k < target_dchisq.size(); k++){
        
        target = target_dchisq[k];
        
        // Initialise box_minmax
        if (k != 0) {
            for (int i = 0; i < np-1; i++) {
                for (int j = 0; j < np-1; j++) {
                    box_minmax[i][j] = {0, 0};
                }
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
            sin_contour[k].push_back(contour[j][0]);
            dm2_contour[k].push_back(contour[j][1]);
        }
        
    }
    
    clock_t endcont = clock();
    clock_t tickscont = endcont - startcont;                    // in n of ticks
    double timecont = tickscont / (double) CLOCKS_PER_SEC;      // make into secs
    
    cout << "Done in " << timecont << "s." << endl;
    
    cout << "Sin contour sizes: " << sin_contour[0].size() << ", " << sin_contour[1].size() << ", " << sin_contour[2].size() << endl;
    cout << "dm2 contour sizes: " << dm2_contour[0].size() << ", " << dm2_contour[1].size() << ", " << dm2_contour[2].size() << endl;
    
    // Now plot it
    TCanvas *canvas = new TCanvas();
    
	TGraph *gr_90 = new TGraph();
    for (int i = 0; i < sin_contour[0].size(); i++) {
        gr_90->SetPoint(i, sin_contour[0][i], dm2_contour[0][i]);
    }
    gr_90->SetMarkerStyle(20);
    gr_90->SetMarkerSize(0.1);
    gr_90->SetMarkerColor(30);
    gr_90->SetLineColor(30);
    
    TGraph *gr_3 = new TGraph();
    for (int i = 0; i < sin_contour[1].size(); i++) {
        gr_3->SetPoint(i, sin_contour[1][i], dm2_contour[1][i]);
    }
    gr_3->SetMarkerStyle(20);
    gr_3->SetMarkerSize(0.1);
    gr_3->SetMarkerColor(38);
    gr_3->SetLineColor(38);
    
    TGraph *gr_5 = new TGraph();
    for (int i = 0; i < sin_contour[2].size(); i++) {
        gr_5->SetPoint(i, sin_contour[2][i], dm2_contour[2][i]);
    }
    gr_5->SetMarkerStyle(20);
    gr_5->SetMarkerSize(0.1);
    gr_5->SetMarkerColor(46);
    gr_5->SetLineColor(46);
    
    TGraph *range = new TGraph();
    range->SetPoint(0, 0.001, 0.01);
    range->SetPoint(1, 1, 100);
    range->SetMarkerColor(0);
    
    gr_90->SetTitle("ICARUS Sensitivity - 30% Uncertainty in Flux; sin^{2}(2#theta); #Delta m^{2} (eV^{2})");
    
    TLegend *legend = new TLegend();
    legend->AddEntry(gr_90, "90% CL", "l");
    legend->AddEntry(gr_3, "3#sigma CL", "l");
    legend->AddEntry(gr_5, "5#sigma CL", "l");
    
    canvas->SetLogy();
    canvas->SetLogx();
    
    gr_90->Draw("AP");
    gr_90->GetXaxis()->SetLimits(0.001, 1);
    gr_90->GetYaxis()->SetRangeUser(0.01, 100);
    canvas->Update();
    
    gr_3->Draw("P same");
    gr_5->Draw("P same");
    legend->Draw();
    range->Draw("P same");
    
    canvas->SaveAs("../Out/LargeError/Sensitivity_LE.root");
    canvas->SaveAs("../Out/LargeError/Sensitivity_LE.png");
    
    
    //// End timer
    clock_t end = clock();
    clock_t ticks = end - start;                        // in n of ticks
    double time = ticks / (double) CLOCKS_PER_SEC;      // make into secs
    
    cout << "Time running: " << time << endl;
    
}