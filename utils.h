#pragma once
#include "RtypesCore.h"
#include "return_TOF_position.h"
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include <string>

using namespace ROOT;

extern int RUN_NUMBER;
extern int BEAM_MOMENTUM;

std::vector<double> extract_SiPM_resolution(
    std::string fname =
        "/home/frantisek/Analysis/configs/T5_reconstruction_values.json");
std::vector<double> convert_resolutions(std::vector<double>);
void write_description();

class Cuts {
  private:
    RVecI T0_ids;
    RVecI T1_ids;
    RVecI T4_ids;
    RVecI HC_ids;
    RVecI T5_ids;

    int T5_board_id;

    double HC_0_threshold;
    double HC_1_threshold;
    double T4_0_threshold;
    double T4_1_threshold;

  public:
    Cuts();

    int get_T5_board();

    RVecI &Get_T5_ids();

    bool hit_T0_T1(const RVecI &bm_time_ids, const RVecI &bm_charge_ids);

    bool did_not_hit_HC(const RVecI &bm_charge_ids, const RVecF &bm_charges);

    bool hit_T4(const RVecI &bm_charge_ids, const RVecF &bm_charges);

    bool hit_T5(const RVecI &mpmt_ids, const RVecI &pmt_ids);
};

class Histograms {
  private:
    std::map<std::string, TH1 *> m_hists_map;
    int m_canvas_default_size_x;
    int m_canvas_default_size_y;
    int m_default_bins;
    double m_default_min;
    double m_default_max;

  public:
    Histograms(int default_can_size_x = 1800, int default_can_size_y = 900,
               int default_bins = 100, double default_min = 0,
               double default_max = 100);
    ~Histograms();

    int get_canvas_default_size(int i);
    void set_canvas_default_size(const int size_x, const int size_y);

    double get_default_range(int i);
    int get_default_bins();
    void set_default_range(const double range_min, const double range_max);
    void set_default_bins(const int bins);

    void book1D(const std::string &name, const std::string &title);
    void book1D(const std::string &name, const std::string &title, int bins,
                double xmin, double xmax);
    void book2D(const std::string &name, const std::string &title);
    void book2D(const std::string &name, const std::string &title, int xbins,
                double xmin, double xmax, int ybins, double ymin, double ymax);
    void fill(const std::string &name, double value);
    void fill(const std::string &name, double x, double y);
    void draw(const std::string &name);
    void draw(const std::string &name, int c_size_x, int c_size_y);
    void print(const std::string &name);
    void print_exclusive(const std::string &name, int c_size_x, int c_size_y);
    void print_all();
    void print_exclusive_log(const std::string &name, int c_size_x,
                             int c_size_y);
    void save_all(const std::string &output_filename);
    TH1 *get_histogram(const std::string &name);
    TH2 *get_histogram_2D(const std::string &name);
    void hist_projectX(const std::string &source_name,
                       const std::string &new_name, int first_y_bin,
                       int last_y_bin);
};
