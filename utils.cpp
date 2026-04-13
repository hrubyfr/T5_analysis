#include "utils.h"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TString.h"
#include "buffer.h"
#include <TLatex.h>
#include <string>

int RUN_NUMBER = 0;
int BEAM_MOMENTUM = 0;

void write_description() {
  TString text = "Run number = " + std::to_string(RUN_NUMBER) +
                 ", Beam momentum = " + std::to_string(BEAM_MOMENTUM) +
                 " MeV/c charged hadron";
  TLatex txt;
  txt.SetTextSize(0.025);
  txt.SetTextAlign(13);
  txt.DrawLatexNDC(0.10, 0.92, text.Data());
}

// Cutting class

Cuts::Cuts()
    : T0_ids{0, 1, 2, 3}, T1_ids{4, 5, 6, 7}, T4_ids{43, 44}, HC_ids{9, 10},
      T5_ids{0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16},

      T5_board_id(132),

      HC_0_threshold(100), HC_1_threshold(100) {}

int Cuts::get_T5_board() { return T5_board_id; }

RVecI &Cuts::Get_T5_ids() { return T5_ids; }

bool Cuts::hit_T0_T1(const RVecI &bm_time_ids, const RVecI &bm_charge_ids) {
  for (auto val : T0_ids) {
    if (!VecOps::Any(bm_time_ids == val)) {
      return false;
    }
  }
  for (auto val : T1_ids) {
    if (!VecOps::Any(bm_time_ids == val)) {
      return false;
    }
  }
  for (auto val : T0_ids) {
    if (!VecOps::Any(bm_charge_ids == val)) {
      return false;
    }
  }
  for (auto val : T1_ids) {
    if (!VecOps::Any(bm_charge_ids == val)) {
      return false;
    }
  }
  return true;
}
bool Cuts::did_not_hit_HC(const RVecI &bm_charge_ids, const RVecF &bm_charges) {
  auto mask_HC0 = (bm_charge_ids == HC_ids[0]);
  auto mask_HC1 = (bm_charge_ids == HC_ids[1]);

  if (VecOps::Any(bm_charges[mask_HC0] > HC_0_threshold) ||
      VecOps::Any(bm_charges[mask_HC1] > HC_1_threshold))
    return false;
  return true;
}
bool Cuts::hit_T4(const RVecI &bm_charge_ids, const RVecF &bm_charges) {
  auto mask_T4_left = (bm_charge_ids == T4_ids[0]);
  auto mask_T4_right = (bm_charge_ids == T4_ids[1]);

  if (VecOps::Any(bm_charges[mask_T4_left] < T4_0_threshold) ||
      VecOps::Any(bm_charges[mask_T4_right] < T4_1_threshold))
    return false;
  return true;
}
bool Cuts::hit_T5(const RVecI &mpmt_ids, const RVecI &pmt_ids) {
  auto mask_board = (mpmt_ids == T5_board_id);
  auto ids_filtered = pmt_ids[mask_board];
  auto half_size = T5_ids.size() / 2;
  for (size_t i = 0; i < half_size; i++) {
    if (VecOps::Any(ids_filtered == T5_ids[i]) &&
        VecOps::Any(ids_filtered == T5_ids[i + half_size]))
      return true;
  }
  return false;
}

// Histograms class

Histograms::Histograms(int default_can_size_x, int default_can_size_y,
                       int default_bins, double default_min, double default_max)
    : m_canvas_default_size_x(default_can_size_x),
      m_canvas_default_size_y(default_can_size_y), m_default_bins(default_bins),
      m_default_min(default_min), m_default_max(default_max) {}

Histograms::~Histograms() {
  for (auto &[name, hist] : m_hists_map) {
    delete hist;
  }
}
int Histograms::get_canvas_default_size(int i) {
  if (i == 0)
    return m_canvas_default_size_x;
  else if (i == 1)
    return m_canvas_default_size_y;
  else
    return -1;
}
void Histograms::set_canvas_default_size(const int size_x, const int size_y) {
  m_canvas_default_size_x = size_x;
  m_canvas_default_size_y = size_y;
}
double Histograms::get_default_range(int i) {
  if (i == 0)
    return m_default_min;
  else if (i == 1)
    return m_default_max;
  else
    return -1;
}
int Histograms::get_default_bins() { return m_default_bins; }
void Histograms::set_default_range(const double min, const double max) {
  m_default_min = min;
  m_default_max = max;
}
void Histograms::set_default_bins(const int bins) { m_default_bins = bins; }

void Histograms::book1D(const std::string &name, const std::string &title) {
  auto hist = new TH1D(name.c_str(), title.c_str(), m_default_bins,
                       m_default_min, m_default_max);
  hist->SetDirectory(nullptr);
  m_hists_map[name] = hist;
}
void Histograms::book1D(const std::string &name, const std::string &title,
                        int bins, double xmin, double xmax) {
  auto hist = new TH1D(name.c_str(), title.c_str(), bins, xmin, xmax);
  hist->SetDirectory(nullptr);
  m_hists_map[name] = hist;
}
void Histograms::book2D(const std::string &name, const std::string &title) {
  auto hist =
      new TH2D(name.c_str(), title.c_str(), m_default_bins, m_default_min,
               m_default_max, m_default_bins, m_default_min, m_default_max);
  hist->SetDirectory(nullptr);
  m_hists_map[name] = hist;
}
void Histograms::book2D(const std::string &name, const std::string &title,
                        int xbins, double xmin, double xmax, int ybins,
                        double ymin, double ymax) {
  auto hist = new TH2D(name.c_str(), title.c_str(), xbins, xmin, xmax, ybins,
                       ymin, ymax);
  hist->SetDirectory(nullptr);
  m_hists_map[name] = hist;
}
void Histograms::fill(const std::string &name, double value) {
  m_hists_map[name]->Fill(value);
}
void Histograms::fill(const std::string &name, double x, double y) {
  m_hists_map[name]->Fill(x, y);
}
void Histograms::draw(const std::string &name) {
  TString c_name = "c_" + name;
  auto canvas = new TCanvas(c_name, c_name, m_canvas_default_size_x,
                            m_canvas_default_size_y);
  canvas->cd();
  m_hists_map[name]->Draw("hist");
}
void Histograms::draw(const std::string &name, int c_size_x, int c_size_y) {
  TString c_name = "c_" + name;
  auto canvas = new TCanvas(c_name, c_name, c_size_x, c_size_y);
  canvas->cd();
  m_hists_map[name]->Draw("hist");
}
void Histograms::print(const std::string &name) {
  TString c_name = "c_" + name;
  auto canvas = new TCanvas(c_name, c_name, m_canvas_default_size_x,
                            m_canvas_default_size_y);
  if (m_hists_map[name] != nullptr) {
    if (m_hists_map[name]->InheritsFrom("TH2")) {
      m_hists_map[name]->SetStats(0);
      m_hists_map[name]->Draw("COLZ");
      write_description();
    } else {
      m_hists_map[name]->SetStats(1);
      m_hists_map[name]->Draw("");
      write_description();
    }
  }
  TString out_name = name + ".png";
  canvas->Print(out_name);
  delete canvas;
}
void Histograms::print_exclusive(const std::string &name, int c_size_x,
                                 int c_size_y) {
  TString c_name = "c_" + name;
  auto canvas = new TCanvas(c_name, c_name, c_size_x, c_size_y);
  canvas->SetRightMargin(0.16);
  if (m_hists_map[name] != nullptr) {
    special_hists.push_back(name);
    if (m_hists_map[name]->InheritsFrom("TH2")) {
      m_hists_map[name]->SetStats(0);
      m_hists_map[name]->Draw("COLZ");
      write_description();
    } else {
      m_hists_map[name]->SetStats(1);
      m_hists_map[name]->Draw("");
      write_description();
    }
  }
  TString out_name = name + ".png";
  canvas->Print(out_name);
  delete canvas;
}

void Histograms::print_all() {
  for (const auto &[name, hist] : m_hists_map) {
    if (std::find(special_hists.begin(), special_hists.end(), name) !=
        special_hists.end())
      continue;
    TString c_name = "c_" + name;
    auto canvas = new TCanvas(c_name, c_name, m_canvas_default_size_x,
                              m_canvas_default_size_y);
    if (hist != nullptr) {
      if (hist->InheritsFrom("TH2")) {
        hist->Draw("COLZ");
        write_description();
      } else {
        hist->Draw("");
        write_description();
      }
    }
    TString out_name = name + ".png";
    canvas->Print(out_name);
    delete canvas;
  }
}
void Histograms::save_all(const std::string &output_name) {
  TString outfile_name = output_name + ".root";
  auto outfile = TFile::Open(outfile_name, "RECREATE");
  for (auto const &[name, hist] : m_hists_map) {
    hist->Write(); // Automatically works for both 1D and 2D
  }
  outfile->Close();
  delete outfile;
}
TH1 *Histograms::get_histogram(const std::string &name) {
  if (m_hists_map.count(name) > 0) {
    return m_hists_map[name];
  }
  std::cerr << "Error: Could not find histogram '" << name << "'\n";
  return nullptr;
}
TH2 *Histograms::get_histogram_2D(const std::string &name) {
  TH1 *hist_ptr = get_histogram(name);
  if (hist_ptr != nullptr) {
    TH2 *hist_2D = dynamic_cast<TH2 *>(hist_ptr);
    return hist_2D;
  }
  std::cerr << "Error: Could not find histogram '" << name
            << "' or is not 2D\n";
  return nullptr;
}
void Histograms::hist_projectX(const std::string &source_name,
                               const std::string &new_name, int first_y_bin,
                               int last_y_bin) {
  // 1. Check if the source histogram actually exists
  if (m_hists_map.count(source_name) == 0) {
    std::cerr << "Error: Source histogram '" << source_name << "' not found!\n";
    return;
  }

  // 2. Safely cast the generic TH1* to a TH2* pointer
  TH2 *h2 = dynamic_cast<TH2 *>(m_hists_map[source_name]);

  // If it was actually a 1D histogram, the cast fails and returns nullptr
  if (h2 == nullptr) {
    std::cerr << "Error: '" << source_name
              << "' is a 1D histogram, cannot ProjectX!\n";
    return;
  }

  // 3. Perform the projection!
  // ROOT automatically creates a new TH1D on the heap and returns the pointer
  TH1D *proj = h2->ProjectionX(new_name.c_str(), first_y_bin, last_y_bin);

  // 4. Detach it from ROOT's memory manager and store it in our dictionary
  proj->SetDirectory(nullptr);
  m_hists_map[new_name] = proj;
}
