#ifndef PROPLOT_H
#define PROPLOT_H

// C++ include 
#include <algorithm>
#include <unordered_map>
#include <string>
#include <iomanip>
// PROfit include 
#include "PROlog.h"
#include "PROconfig.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROMCMC.h"
#include "PROtocall.h"
#include "PROseed.h"
#include "PROcess.h"
#include "PROversion.h"

// Root includes
#include "TAttLine.h"
#include "TAttMarker.h"
#include "THStack.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRatioPlot.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TLine.h"
namespace PROfit{

    class ROOTFileWriter
    {
    public:
        ROOTFileWriter(const ROOTFileWriter&) = delete;
        ROOTFileWriter& operator=(const ROOTFileWriter&) = delete;

        // Get the singleton instance
        static ROOTFileWriter& instance()
        {
            static ROOTFileWriter instance;
            return instance;
        }

        // Open a ROOT file
        void open(const std::string & filename, const char * option = "RECREATE")
        {
            // Close any previously open file before opening a new one
            if(file_ && file_->IsOpen())
            {
                log<LOG_WARNING>(L"ROOTFileWriter: File already open. Closing previous file.");
                close();
            }
            
            // Persist the current directory to restore it later
            TDirectory * current_dir = gDirectory;

            // Open the new file
            file_ = new TFile(filename.c_str(), option);

            // Check if the file was opened successfully
            if(!file_ || file_->IsZombie())
            {
                log<LOG_ERROR>(L"ROOTFileWriter: Failed to open ROOT file: %1%") % filename.c_str();
                if (file_)
                {
                    delete file_;
                    file_ = nullptr;
                }
                if(current_dir) current_dir->cd();
                throw std::runtime_error("Failed to open ROOT file: " + filename);
            }

            // Set the file name for reference
            filename_ = filename;
            log<LOG_INFO>(L"ROOTFileWriter: Opened ROOT file: %1%") % filename.c_str();

            // Restore the original directory
            if(current_dir) current_dir->cd();

            // Create a TTree to store 1D histogram data
            tree_hist1d_ = new TTree("hist1d", "1D histogram data");
            tree_hist1d_->SetDirectory(nullptr);
            tree_hist1d_->Branch("variable", &hist1d_variable_, "variable/i");
            tree_hist1d_->Branch("mode", &hist1d_mode_, "mode/i");
            tree_hist1d_->Branch("detector", &hist1d_detector_, "detector/i");
            tree_hist1d_->Branch("channel", &hist1d_channel_, "channel/i");
            tree_hist1d_->Branch("subchannel", &hist1d_subchannel_, "subchannel/i");
            tree_hist1d_->Branch("prefix", &hist1d_prefix_);
            tree_hist1d_->Branch("bin_index", &hist1d_bin_index_, "bin_index/I");
            tree_hist1d_->Branch("bin_center", &hist1d_bin_center_, "bin_center/D");
            tree_hist1d_->Branch("bin_low_edge", &hist1d_bin_low_edge_, "bin_low_edge/D");
            tree_hist1d_->Branch("bin_high_edge", &hist1d_bin_high_edge_, "bin_high_edge/D");
            tree_hist1d_->Branch("bin_content", &hist1d_bin_content_, "bin_content/D");
            tree_hist1d_->Branch("bin_error",   &hist1d_bin_error_,   "bin_error/D");

            // Create a TTree to store error band data
            tree_errorband_ = new TTree("errorband", "Error band data");
            tree_errorband_->SetDirectory(nullptr);
            tree_errorband_->Branch("variable", &errband_variable_, "variable/i");
            tree_errorband_->Branch("mode", &errband_mode_, "mode/i");
            tree_errorband_->Branch("detector", &errband_detector_, "detector/i");
            tree_errorband_->Branch("channel", &errband_channel_, "channel/i");
            tree_errorband_->Branch("subchannel", &errband_subchannel_, "subchannel/i");
            tree_errorband_->Branch("prefix", &errband_prefix_);
            tree_errorband_->Branch("point_index", &errband_point_index_, "point_index/I");
            tree_errorband_->Branch("x_value", &errband_x_value_, "x_value/D");
            tree_errorband_->Branch("y_value", &errband_y_value_, "y_value/D");
            tree_errorband_->Branch("error_y_low", &errband_error_y_low_, "error_y_low/D");
            tree_errorband_->Branch("error_y_high", &errband_error_y_high_, "error_y_high/D");

            // Create a TTree to store fractional systematic data
            tree_frac_syst_ = new TTree("frac_syst", "Fractional systematic data");
            tree_frac_syst_->SetDirectory(nullptr);
            tree_frac_syst_->Branch("mode", &frac_syst_mode_, "mode/i");
            tree_frac_syst_->Branch("detector", &frac_syst_detector_, "detector/i");
            tree_frac_syst_->Branch("channel", &frac_syst_channel_, "channel/i");
            tree_frac_syst_->Branch("tag", &frac_syst_tag_);
            tree_frac_syst_->Branch("systname", &frac_syst_systname_);
            tree_frac_syst_->Branch("bin_index", &frac_syst_bin_index_, "bin_index/I");
            tree_frac_syst_->Branch("bin_center", &frac_syst_bin_center_, "bin_center/D");
            tree_frac_syst_->Branch("bin_low_edge", &frac_syst_bin_low_edge_, "bin_low_edge/D");
            tree_frac_syst_->Branch("bin_high_edge", &frac_syst_bin_high_edge_, "bin_high_edge/D");
            tree_frac_syst_->Branch("bin_content", &frac_syst_bin_content_, "bin_content/D");
        }

        // Close the current file
        void close()
        {
            if(file_ && file_->IsOpen())
            {
                // Persist the current directory to restore it later
                TDirectory * current_dir = gDirectory;

                file_->cd();

                if(tree_hist1d_)
                {
                    tree_hist1d_->SetDirectory(file_); 
                    tree_hist1d_->Write("hist1d", TObject::kOverwrite);
                    tree_hist1d_ = nullptr;
                }
                if(tree_errorband_)
                {
                    tree_errorband_->SetDirectory(file_);
                    tree_errorband_->Write("errorband", TObject::kOverwrite);
                    tree_errorband_ = nullptr;
                }
                if(tree_frac_syst_)
                {
                    tree_frac_syst_->SetDirectory(file_);
                    tree_frac_syst_->Write("frac_syst", TObject::kOverwrite);
                    tree_frac_syst_ = nullptr;
                }

                file_->Close();
                //delete file_;
                file_ = nullptr;
            }
            filename_.clear();
        }

        // Check if the file is already open
        bool is_open() const
        {
            return file_ && file_->IsOpen();
        }

        // Set the variable number
        void set_variable(UInt_t variable)
        {
            hist1d_variable_ = variable;
            errband_variable_ = variable;
        }

        // Static member function to get the variable number from the file name
        static int extract_variable_number(const std::string& filename)
        {
            // Look for pattern "Variable_XX"
            size_t pos = filename.find("Variable_");
            if(pos == std::string::npos) return -1;
            
            pos += 9;
            
            // Extract digits
            size_t end = pos;
            while (end < filename.length() && std::isdigit(filename[end]))
                end++;
            
            if (end > pos)
                return std::stoi(filename.substr(pos, end - pos));
            
            return -1;
        }

        // Populate a row in the hist1d TTree
        void fill_hist1d(const TH1D& hist, const std::string& prefix, 
                            size_t mode, size_t det, size_t channel, size_t subchannel)
        {
            if(!is_open() || !tree_hist1d_) return;

            for(int i = 1; i <= hist.GetNbinsX(); ++i)
            {
                hist1d_mode_ = mode;
                hist1d_detector_ = det;
                hist1d_channel_ = channel;
                hist1d_subchannel_ = subchannel;
                *hist1d_prefix_ = prefix;
                hist1d_bin_index_ = i - 1;
                hist1d_bin_center_ = hist.GetBinCenter(i);
                hist1d_bin_low_edge_ = hist.GetBinLowEdge(i);
                hist1d_bin_high_edge_ = hist.GetBinLowEdge(i) + hist.GetBinWidth(i);
                hist1d_bin_content_ = hist.GetBinContent(i);
                hist1d_bin_error_   = hist.GetBinError(i);
                if(tree_hist1d_) tree_hist1d_->Fill();
            }
        }

        // Populate rows in the error band TTree
        void fill_errorband(const TGraphAsymmErrors* errband, const std::string& prefix,
                            size_t mode, size_t det, size_t channel, size_t subchannel)
        {
            if(!errband || !is_open() || !tree_errorband_) return;

            for(int i = 0; i < errband->GetN(); ++i)
            {
                double x, y;
                errband->GetPoint(i, x, y);
                errband_mode_ = mode;
                errband_detector_ = det;
                errband_channel_ = channel;
                errband_subchannel_ = subchannel;
                *errband_prefix_ = prefix;
                errband_point_index_ = i;
                errband_x_value_ = x;
                errband_y_value_ = y;
                errband_error_y_low_ = errband->GetErrorYlow(i);
                errband_error_y_high_ = errband->GetErrorYhigh(i);
                if(tree_errorband_) tree_errorband_->Fill();
            }
        }

        // Populate rows in the fractional systematic TTree
        void fill_fractional_systematic(const TH1F* hist, const std::string& tag,
                                const std::string& systname,
                                size_t mode, size_t det, size_t channel)
        {
            if(!hist || !is_open() || !tree_frac_syst_) return;
            
            for(int i = 1; i <= hist->GetNbinsX(); ++i)
            {
                frac_syst_mode_ = mode;
                frac_syst_detector_ = det;
                frac_syst_channel_ = channel;
                *frac_syst_tag_ = tag;
                *frac_syst_systname_ = systname;
                frac_syst_bin_index_ = i - 1;
                frac_syst_bin_center_ = hist->GetBinCenter(i);
                frac_syst_bin_low_edge_ = hist->GetBinLowEdge(i);
                frac_syst_bin_high_edge_ = hist->GetBinLowEdge(i) + hist->GetBinWidth(i);
                frac_syst_bin_content_ = hist->GetBinContent(i);
                if(tree_frac_syst_) tree_frac_syst_->Fill();
            }
        }

    private:
        ROOTFileWriter() : file_(nullptr), 
            tree_hist1d_(nullptr), hist1d_prefix_(new std::string()),
            tree_errorband_(nullptr), errband_prefix_(new std::string()),
            tree_frac_syst_(nullptr), frac_syst_tag_(new std::string()), frac_syst_systname_(new std::string()) {}

        TFile* file_;
        std::string filename_;

        // Variables for hist1d TTree
        TTree* tree_hist1d_;
        UInt_t hist1d_variable_;
        UInt_t hist1d_mode_;
        UInt_t hist1d_detector_;
        UInt_t hist1d_channel_;
        UInt_t hist1d_subchannel_;
        std::string * hist1d_prefix_;
        Int_t hist1d_bin_index_;
        Double_t hist1d_bin_center_;
        Double_t hist1d_bin_low_edge_;
        Double_t hist1d_bin_high_edge_;
        Double_t hist1d_bin_content_;
        Double_t hist1d_bin_error_;

        // Variables for error band TTree
        TTree* tree_errorband_;
        UInt_t errband_variable_;
        UInt_t errband_mode_;
        UInt_t errband_detector_;
        UInt_t errband_channel_;
        UInt_t errband_subchannel_;
        std::string * errband_prefix_;
        Int_t errband_point_index_;
        Double_t errband_x_value_;
        Double_t errband_y_value_;
        Double_t errband_error_y_low_;
        Double_t errband_error_y_high_;

        // Variables for fractional systematic TTree
        TTree* tree_frac_syst_;
        UInt_t frac_syst_mode_;
        UInt_t frac_syst_detector_;
        UInt_t frac_syst_channel_;
        std::string * frac_syst_tag_;
        std::string * frac_syst_systname_;
        Int_t frac_syst_bin_index_;
        Double_t frac_syst_bin_center_;
        Double_t frac_syst_bin_low_edge_;
        Double_t frac_syst_bin_high_edge_;
        Double_t frac_syst_bin_content_;

        ~ROOTFileWriter()
        { 
            // String pointers deleted after file is closed
            delete hist1d_prefix_;
            delete errband_prefix_;
            delete frac_syst_tag_;
            delete frac_syst_systname_;
        }
    };

    struct PlotBounds {
        float xmin = -9999;
        float xmax = -9999;
        float ymin = -9999;
        float ymax = -9999;
        float ratmin = -9999;
        float ratmax = -9999;

        int Load(std::map<std::string, float> bound_list){
            log<LOG_INFO>(L"%1% || Loading Bounds for plot_channels ") % __func__;
            for(const auto &[bound_name, value]: bound_list) {
                log<LOG_INFO>(L"%1% || --on bound %2% val %3% ") % __func__ % bound_name.c_str() % value;
                if(bound_name == "xmin") {
                    xmin=value;
                }else if(bound_name == "xmax") {
                    xmax=value;
                }else if(bound_name == "ymin") {
                    ymin=value;
                }else if(bound_name == "ymax") {
                    ymax=value;
                }else if(bound_name == "ratmin") {
                    ratmin=value;
                }else if(bound_name == "ratmax") {
                    ratmax=value;
                }else{
                    log<LOG_ERROR>(L"%1% || ERROR! you passed a plot-bounds string that is not allowed (%2%). Needs to be xmin,xmax,ymin,ymax,ratmin,ratmax. ") % __func__ % bound_name.c_str();
                    throw std::invalid_argument(std::string("Invalid plot-bounds ")+bound_name );
                }
            }
        return 1;
        };
        bool hasBound(std::string bound_name){
                if(bound_name == "xmin") {
                    return xmin!=-9999? true : false;
                }else if(bound_name == "xmax") {
                    return xmax!=-9999? true : false;
                }else if(bound_name == "ymin") {
                    return ymin!=-9999? true : false;
                }else if(bound_name == "ymax") {
                    return ymax!=-9999? true : false;
                }else if(bound_name == "ratmin") {
                    return ratmin!=-9999? true : false;
                }else if(bound_name == "ratmax") {
                    return ratmax!=-9999? true : false;
                }else{
                    log<LOG_ERROR>(L"%1% || ERROR! you passed a plot-bounds string that is not allowed (%2%). Needs to be ymax,ratmin,ratmax. ") % __func__ % bound_name.c_str();
                    throw std::invalid_argument(std::string("Invalid plot-bounds ")+bound_name );
                }
        return false;
        };
        float getBound(std::string bound_name){
                if(bound_name == "xmin") {
                    return xmin;
                }else if(bound_name == "xmax") {
                    return xmax;
                }else if(bound_name == "ymin") {
                    return ymin;
                }else if(bound_name == "ymax") {
                    return ymax;
                }else if(bound_name == "ratmin") {
                    return ratmin;
                }else if(bound_name == "ratmax") {
                    return ratmax;
                }else{
                    log<LOG_ERROR>(L"%1% || ERROR! you passed a plot-bounds string that is not allowed (%2%). Needs to be ymax,ratmin,ratmax. ") % __func__ % bound_name.c_str();
                    throw std::invalid_argument(std::string("Invalid plot-bounds ")+bound_name );
                }
        return -999;
        };
    };
    
    enum class PlotOptions {
        Default = 0,
        CVasStack = 1 << 0,
        AreaNormalized = 1 << 1,
        BinWidthScaled = 1 << 2,
        DataMCRatio = 1 << 3,
        DataPostfitRatio = 1 << 4,
    };

    inline PlotOptions operator|(PlotOptions a, PlotOptions b) {
        return static_cast<PlotOptions>(static_cast<int>(a) | static_cast<int>(b));
    }

    inline PlotOptions operator|=(PlotOptions &a, PlotOptions b) {
        return a = a | b;
    }

    inline PlotOptions operator&(PlotOptions a, PlotOptions b) {
        return static_cast<PlotOptions>(static_cast<int>(a) & static_cast<int>(b));
    }

    inline PlotOptions operator&=(PlotOptions &a, PlotOptions b) {
        return a = a & b;
    }


    void plot_detector_ratios(const PROconfig &config, std::vector<TH1D> data_hists,  std::vector<TH1D> cv_hists, std::optional<PROerrorbar> errband, std::vector<TH1D> bf_hists, std::optional<PROerrorbar> posterrband, TH2D &pre_corr, TH2D &post_corr, std::string filename, int var_index = 0);

    void plot_channels(const std::string &filename, const PROconfig &config, std::optional<PROspec> cv, std::optional<PROspec> best_fit, std::optional<PROdata> data, std::optional<PROerrorbar> errband, std::optional<PROerrorbar> posterrband, std::optional<PROsyst> pre_allcovsyst, std::optional<PROsyst> post_allcovsyst, std::vector<TPaveText> &texts, PlotBounds &bounds, PlotOptions opt = PlotOptions::Default, int var_index = 0);

    //some helper functions for PROplot
    std::map<std::string, std::unique_ptr<TH1D>> getCV1DHists(const PROspec & spec, const PROconfig& inconfig, bool scale = false, int var_index = 0);
    std::map<std::string, std::unique_ptr<TH2D>> getCV2DHists(const PROspec & spec, const PROconfig& inconfig, bool scale = false, int var_index = 0);
    std::map<std::string, std::unique_ptr<TH2D>> covarianceTH2D(const PROsyst &syst, const PROconfig &config, const PROspec &cv);
    std::map<std::string, std::vector<std::pair<std::unique_ptr<TGraph>,std::unique_ptr<TGraph>>>> getSplineGraphs(const PROsyst &systs, const PROconfig &config);
    PROerrorbar getErrorBand(const PROconfig &config, const PROpeller &prop, const PROsyst &syst, const PROmodel &model, const PROspec &cv_spec, const Eigen::VectorXf &cvparams,bool scale=false, int other_index=0);

    int plotPriorFractionalSystematicBreakdown(const PROconfig &config, const PROspec &spec, const PROsyst &allsplinesyst, std::string filename, int var_index = 0);

    template<class T, class P>
        PROerrorbar getMCMCErrorBand(Metropolis<T, P> met, size_t burnin, size_t iterations, const PROconfig &config, const PROpeller &prop, PROmetric &metric, const Eigen::VectorXf &best_fit, std::vector<TH1D> &posteriors, Eigen::MatrixXf &post_covar,  bool scale = false,int var_index=0) {
            for(size_t i = 0; i < metric.GetSysts().GetNSplines(); ++i)
                posteriors.emplace_back("", (";"+config.m_mcgen_variation_plotname_map.at(metric.GetSysts().spline_names[i])).c_str(), 60, -3, 3);

            Eigen::VectorXf cv = FillSpectra(config, prop, metric.GetSysts(), metric.GetModel(), best_fit, true, var_index).Spec();
            Eigen::MatrixXf L; 
            if(metric.GetSysts().GetNCovar() > 0) L = metric.GetSysts().DecomposeFractionalCovariance(config, cv);
            else L = Eigen::MatrixXf::Zero(config.m_num_variable_bins_total_collapsed[var_index], config.m_num_variable_bins_total_collapsed[var_index]);
            std::normal_distribution<float> nd;
            Eigen::VectorXf throws = Eigen::VectorXf::Constant(config.m_num_variable_bins_total_collapsed[var_index], 0);

            int nspline = metric.GetSysts().GetNSplines();
            int nphys = metric.GetModel().nparams;
            Eigen::VectorXf splines_bf = best_fit.segment(nphys, nspline);
            post_covar = Eigen::MatrixXf::Constant(nspline, nspline, 0);
            size_t accepted = 0;
            std::vector<Eigen::VectorXf> specs;
            const auto action = [&](const Eigen::VectorXf &value) {
                accepted += 1;
                for(size_t i = 0; i < config.m_num_variable_bins_total_collapsed[var_index]; ++i)
                    throws(i) = nd(PROseed::global_rng);
                specs.push_back(CollapseMatrix(config, FillSpectra(config, prop, metric.GetSysts(), metric.GetModel(), value, true,var_index).Spec())+L*throws);
                for(size_t i = 0; i < metric.GetSysts().GetNSplines(); ++i)
                    posteriors[i].Fill(value(i+nphys));
                Eigen::VectorXf splines = value.segment(nphys, nspline);
                Eigen::VectorXf diff = splines-splines_bf;
                post_covar += diff * diff.transpose();
                bool print_autocorrelation_values = false;
                if(print_autocorrelation_values){
                    log<LOG_INFO>(L"%1% || AUTO  %2% : %3%") % __func__ % accepted % value;
                }
            };
            met.run(burnin, iterations, action);
            post_covar /= accepted;
            log<LOG_INFO>(L"%1% || Acceptance rate %2%") % __func__ % ((float)accepted / iterations);

            cv = CollapseMatrix(config, cv);

            std::vector<float> centers;
            size_t global_channel_index = 0;
            for(size_t mode = 0; mode < config.m_num_modes; ++mode) {
                for(size_t det = 0; det < config.m_num_detectors; ++det) {
                    for(size_t channel = 0; channel < config.m_num_channels; ++channel) {
                        std::vector<float> tedges =  config.GetChannelVariableBins(global_channel_index, var_index).Edges();
                        global_channel_index++;
                        for(size_t p=0; p<tedges.size(); p++){
                            if(p<tedges.size()-1){
                                centers.push_back((tedges[p+1]+tedges[p])/2.0);
                            }
                        }

                    }
                }
            }

            PROerrorbar ebar(cv.size());
            for(int i = 0; i < cv.size(); ++i) {
                std::vector<float> binconts(specs.size());
                for(size_t j = 0; j < specs.size(); ++j) {
                    binconts[j] = specs[j](i);
                }
                float scale_factor = scale ? 1.0/config.collapsed_bin_widths.at(var_index)(i) :  1.0;
                if(std::isnan(scale_factor)) scale_factor = 1;
                std::sort(binconts.begin(), binconts.end());
                float ehi = std::abs((binconts[int(0.840*specs.size())] - cv(i))*scale_factor);
                float elo = std::abs((cv(i) - binconts[int(0.160*specs.size())])*scale_factor);
                ebar.error_up(i) =  ehi;
                ebar.error_down(i) =  elo;
                ebar.error_point(i) = cv(i)*scale_factor;
                log<LOG_DEBUG>(L"%1% || ErrorBand bin %2% %3% %4% %5% %6% ") % __func__ % i % cv(i) % ehi % elo % scale_factor ;
            }
            return ebar;
        }


};

#endif
