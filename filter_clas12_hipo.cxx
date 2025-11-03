/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  ROOT Macro to filter clas 12 data from hipo files
///
///  Stefan Diehl  (sdiehl@jlab.org)
///
///  export CLAS12ROOT="/home/sdiehl/clas12root"            
///  export PATH="$PATH":"$CLAS12ROOT/bin"
///  export RCDB_HOME=${CLAS12ROOT}/rcdb
///  export CCDB_HOME=${CLAS12ROOT}/ccdb
///  source ${CCDB_HOME}/environment.bash
///  export QADB=${CLAS12ROOT}/clasqaDB
///  clas12root
///
///  To execute run clas12root and compile macro with .L filter_clas12_hipo.cxx++
///
///  Then call function:  filter_clas12_hipo("nSidis_005043.hipo", "inb_005043_test.root", 5043)
///                       filter_clas12_hipo2("nSidis_005487.hipo", "outb_005487_test.root", 5487)
//
///  FD eid cuts:  0 = negative charge,  1 = PID default, 2 = PID + charge + EC inner vs outer,  3 = PID + charge + EC samp. frac.,  4 = PID + charge + PCAL hit pos. fiducial,
///                5 = PID + charge + DC fiducial region 1,  6 = PID + charge + DC fiducial region 3,  7 = PID + charge + DC z vertex,
///                8 = all cuts passed, except the one based on the shown histogram, 9 = all cuts,  10 = inverse PID cut
///
///  FD hadron ID cuts:  0 = PID default,  1 = charge,   2 = DC fiducial region 1,    3 = DC fiducial region 3,   8 = delta vz, 9 = all
///                      --> cuts 3 to 8 contain charge and DC fiducial cut as basis, for neutrons, the DC fiducial cut is not applied
///                      --> For negative pions and Kaons also a a cut, requireing, that the particle is not a correctly identified electron is applied to cut 1 - 7
///                      --> For negative pion also a EC inner vs outer cut is applied to reject electrons to cut 1 - 7
///
///                      Proton = + 0,  Neutron = + 10,  Pip = + 20,  Pim = + 30,  Kp = + 40,  Km = + 50
///
///  CD hadron ID:  0 = PID default,   1 = charge,   2 = PID + charge + beta vs p,   3 = PID + charge + maximum probability,   4 = PID + charge + z vertex,  5 = all cuts
///
///  Photon ID:  0 = PID default,  1 = charge,  2 = beta cut,  3 = PCAL hit position fiducial cut,  4 = all cuts
///
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <TChain.h>
#include <TCutG.h>
#include <fstream>
#include <bitset>
using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "clas12reader.h"
#include "jsonFileMerger.h"

using namespace clas12;

using namespace ROOT::Math;

// #include "dataset/weight3FinValFTOF/TMVAClassification_DNN.class.C"

// #include "TMVA/Reader.h"

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

int process_Events = 10000;
//int process_Events = 100000;

// float Ebeam = 6.42313;
// float Ebeam = 7.54626;
// float Ebeam = 2.22193;

 float Ebeam = 10.6041;
// float Ebeam = 10.1998;

/// select the field polarity (relevant for fiducial cuts, sampling fraction and vertex cuts)
/// set only one as true!

bool inbending = true;
bool outbending = false;
bool spring2019 = false;

///

bool simulation = false;


///

bool use_qa_cuts = false;

///

bool use_FT = false;
bool use_FD = true;
bool use_CD = true;

/////////////////////////////////////////////////////////////////////////////////
/// fill output variables (beta, LTCC, ML, RICH, PCAL and or DC coordinates:

bool writebeta  = false;
bool writeML    = false;
bool writeLTCC  = false;
bool userich    = true;
bool fill_coord_ele  = false;
bool fill_coord_prot = true;
bool fill_coord_pip  = false;
bool fill_coord_pim  = false;
bool fill_coord_phot = false;
bool fill_edge_ele  = true;
bool fill_edge_prot = false;
bool fill_edge_pip  = false;
bool fill_edge_pim  = false;


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////
/// select the type of particle ID:

/// FD:

bool use_own_PID_electron = true; // true
bool use_own_PID_photon = true;   // true

/// evb PID is used as basis for p < 2 GeV, select if all neutrals should be considered for p > 2 GeV:

bool cut_extend_neutron = true;

/// select onyl 1 of the following (if nothing is selected then use evb PID):

bool use_fiducial_charged = true; // true     // only fiducial and dvz cuts (automatically inlcuded in the others)
bool cut_maximum_probability_charged = false;
bool population_weighting = false;
bool cut_beta_vs_p_charged = false;
bool cut_deltabeta_charged = false;
bool cut_tofmass_charged = false;

////////////////////////////////////////

// CD:

/// select only 1 of the following (if nothing is selected then use evb PID):

bool use_CD_fiducial_charged = true; // only fiducial and dvz cuts (automatically inlcuded in the others)
bool CD_cut_maximum_probability_charged = false;
bool population_weighting_CD = false;
bool CD_cut_beta_vs_p_charged = false;

/// //////////////////////////////////////////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////////////////////////////////////////


// use default, not corrected beta?

bool show_FD_eid_statistics = true;
bool show_charged_ID_statistics = true;
bool show_photon_ID_statistics = true;

bool fill_electron_pid_histograms = true;
bool fill_hadron_pid_histograms = true;
bool fill_photon_pid_histograms = true;

const static int BUFFER = 60; // increased from 40

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// common variables
///

Float_t Pival = 3.14159265359;

Float_t m_e = 0.000511;
Float_t m_p = 0.93827;
Float_t m_n = 0.9396;
Float_t m_pip = 0.1396;
Float_t m_pim = 0.1396;
Float_t m_Kp = 0.4937;
Float_t m_Km = 0.4937;
Float_t c = 299792458;

/// /////////////////////////////////////////////////////////////////////////////////////
///  output file:

TFile *out;

char name[200];
char title[200];

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// general particle properties and additional particle bank properties for non identified particles:

short TYPE, Helic, Helic_raw;
int NRUN, NEVENT;
long long int TRG;
float EVNTime, BCG, STTime, RFTime;

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// general particle properties and additional particle bank properties for non identified particles:

int Npart;

float part_px[BUFFER], part_py[BUFFER], part_pz[BUFFER];
float part_beta[BUFFER], part_chi2pid[BUFFER];
float part_p[BUFFER], part_theta[BUFFER], part_phi[BUFFER];
float part_vx[BUFFER], part_vy[BUFFER], part_vz[BUFFER];
short part_charge[BUFFER], part_status[BUFFER];
int part_pid[BUFFER];

int ele_detect[BUFFER], prot_detect[BUFFER], neutr_detect[BUFFER], pip_detect[BUFFER], pim_detect[BUFFER], Kp_detect[BUFFER], Km_detect[BUFFER], phot_detect[BUFFER];
int ele_sector[BUFFER], prot_sector[BUFFER], neutr_sector[BUFFER], pip_sector[BUFFER], pim_sector[BUFFER], Kp_sector[BUFFER], Km_sector[BUFFER], phot_sector[BUFFER];
float ele_chi2pid[BUFFER], prot_chi2pid[BUFFER], neutr_chi2pid[BUFFER], pip_chi2pid[BUFFER], pim_chi2pid[BUFFER], Kp_chi2pid[BUFFER], Km_chi2pid[BUFFER];
float ele_dcx1[BUFFER], pip_dcx1[BUFFER], pim_dcx1[BUFFER];
float ele_dcy1[BUFFER], pip_dcy1[BUFFER], pim_dcy1[BUFFER];
float ele_dcx2[BUFFER], pip_dcx2[BUFFER], pim_dcx2[BUFFER];
float ele_dcy2[BUFFER], pip_dcy2[BUFFER], pim_dcy2[BUFFER];
float ele_dcx3[BUFFER], pip_dcx3[BUFFER], pim_dcx3[BUFFER];
float ele_dcy3[BUFFER], pip_dcy3[BUFFER], pim_dcy3[BUFFER];
float prot_dcx1[BUFFER], prot_dcy1[BUFFER], prot_dcz1[BUFFER];
float ele_pcalv[BUFFER], phot_pcalv[BUFFER], ele_pcalw[BUFFER], phot_pcalw[BUFFER];

float ele_dcedge1[BUFFER], prot_dcedge1[BUFFER], pip_dcedge1[BUFFER], pim_dcedge1[BUFFER];
float ele_dcedge2[BUFFER], prot_dcedge2[BUFFER], pip_dcedge2[BUFFER], pim_dcedge2[BUFFER];
float ele_dcedge3[BUFFER], prot_dcedge3[BUFFER], pip_dcedge3[BUFFER], pim_dcedge3[BUFFER];


/// ///////////////////////////////
/// MC:

float MC_Ebeam, MC_weight;
int MC_helicity, MC_Npart;

int partMC_pid[BUFFER];
float partMC_px[BUFFER], partMC_py[BUFFER], partMC_pz[BUFFER];
float partMC_p[BUFFER], partMC_theta[BUFFER], partMC_phi[BUFFER];
float partMC_vx[BUFFER], partMC_vy[BUFFER], partMC_vz[BUFFER];
int partMC_mother[BUFFER];

int partLUND_pid[BUFFER], partLUND_mother[BUFFER];
float partLUND_mass[BUFFER], partLUND_E[BUFFER];
float partLUND_px[BUFFER], partLUND_py[BUFFER], partLUND_pz[BUFFER];
float partLUND_p[BUFFER], partLUND_theta[BUFFER], partLUND_phi[BUFFER];
float partLUND_vx[BUFFER], partLUND_vy[BUFFER], partLUND_vz[BUFFER];

/// //////////////////////////////

int part_Cal_PCAL_sector[BUFFER], part_Cal_ECin_sector[BUFFER], part_Cal_ECout_sector[BUFFER];
float part_Cal_PCAL_energy[BUFFER], part_Cal_ECin_energy[BUFFER], part_Cal_ECout_energy[BUFFER], part_Cal_energy_total[BUFFER];
float part_Cal_PCAL_time[BUFFER], part_Cal_ECin_time[BUFFER], part_Cal_ECout_time[BUFFER];
float part_Cal_PCAL_path[BUFFER], part_Cal_ECin_path[BUFFER], part_Cal_ECout_path[BUFFER];
float part_Cal_PCAL_x[BUFFER], part_Cal_PCAL_y[BUFFER], part_Cal_PCAL_z[BUFFER];
float part_Cal_ECin_x[BUFFER], part_Cal_ECin_y[BUFFER], part_Cal_ECin_z[BUFFER];
float part_Cal_ECout_x[BUFFER], part_Cal_ECout_y[BUFFER], part_Cal_ECout_z[BUFFER];
float part_Cal_PCAL_lu[BUFFER], part_Cal_PCAL_lv[BUFFER], part_Cal_PCAL_lw[BUFFER];
float part_Cal_ECin_lu[BUFFER], part_Cal_ECin_lv[BUFFER], part_Cal_ECin_lw[BUFFER];
float part_Cal_ECout_lu[BUFFER], part_Cal_ECout_lv[BUFFER], part_Cal_ECout_lw[BUFFER];

float part_Cal_PCAL_m2u[BUFFER], part_Cal_ECin_m2u[BUFFER], part_Cal_ECout_m2u[BUFFER];
float part_Cal_PCAL_m2v[BUFFER], part_Cal_ECin_m2v[BUFFER], part_Cal_ECout_m2v[BUFFER];
float part_Cal_PCAL_m2w[BUFFER], part_Cal_ECin_m2w[BUFFER], part_Cal_ECout_m2w[BUFFER];
float part_Cal_PCAL_m3u[BUFFER], part_Cal_ECin_m3u[BUFFER], part_Cal_ECout_m3u[BUFFER];
float part_Cal_PCAL_m3v[BUFFER], part_Cal_ECin_m3v[BUFFER], part_Cal_ECout_m3v[BUFFER];
float part_Cal_PCAL_m3w[BUFFER], part_Cal_ECin_m3w[BUFFER], part_Cal_ECout_m3w[BUFFER];

int part_CC_HTCC_sector[BUFFER], part_CC_HTCC_nphe[BUFFER];
float part_CC_HTCC_time[BUFFER], part_CC_HTCC_path[BUFFER];
float part_CC_HTCC_theta[BUFFER], part_CC_HTCC_phi[BUFFER];

int part_CC_LTCC_sector[BUFFER], part_CC_LTCC_nphe[BUFFER];
float part_CC_LTCC_time[BUFFER], part_CC_LTCC_path[BUFFER];
float part_CC_LTCC_theta[BUFFER], part_CC_LTCC_phi[BUFFER];

int part_FTOF_layer[BUFFER];
int part_FTOF_sector_layer1[BUFFER], part_FTOF_sector_layer2[BUFFER], part_FTOF_sector_layer3[BUFFER];
int part_FTOF_component_layer1[BUFFER], part_FTOF_component_layer2[BUFFER], part_FTOF_component_layer3[BUFFER];
float part_FTOF_energy[BUFFER], part_FTOF_time[BUFFER], part_FTOF_path[BUFFER];
float part_FTOF_energy_layer1[BUFFER], part_FTOF_time_layer1[BUFFER], part_FTOF_path_layer1[BUFFER];
float part_FTOF_energy_layer3[BUFFER], part_FTOF_time_layer3[BUFFER], part_FTOF_path_layer3[BUFFER];

int part_CTOF_component[BUFFER];
float part_CTOF_energy[BUFFER], part_CTOF_time[BUFFER], part_CTOF_path[BUFFER];

int part_CND_component[BUFFER];
float part_CND_energy[BUFFER], part_CND_time[BUFFER], part_CND_path[BUFFER];

float part_FT_energy[BUFFER], part_FT_radius[BUFFER], part_FTHODO_time[BUFFER], part_FTHODO_path[BUFFER], part_FTCAL_time[BUFFER], part_FTCAL_path[BUFFER];
float part_FTCAL_x[BUFFER], part_FTCAL_y[BUFFER], part_FTCAL_z[BUFFER];
float part_FTTRK_x[BUFFER], part_FTTRK_y[BUFFER], part_FTTRK_z[BUFFER];
float part_FTHODO_x[BUFFER], part_FTHODO_y[BUFFER], part_FTHODO_z[BUFFER];

int part_DC_sector[BUFFER], part_DC_index[BUFFER];
int part_DC_Track_NDF[BUFFER], part_DC_Track_status[BUFFER];
float part_DC_Track_chi2[BUFFER];
float part_DC_c1x[BUFFER], part_DC_c1y[BUFFER], part_DC_c1z[BUFFER];
float part_DC_c2x[BUFFER], part_DC_c2y[BUFFER], part_DC_c2z[BUFFER];
float part_DC_c3x[BUFFER], part_DC_c3y[BUFFER], part_DC_c3z[BUFFER];

float part_DC_edge1[BUFFER], part_DC_edge2[BUFFER], part_DC_edge3[BUFFER];

int part_RICH_best[BUFFER];

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// partcile variables for identified particles:

int e_count, p_count, n_count, pip_count, pim_count, Kp_count, Km_count, g_count;
int e_MCcount, p_MCcount, n_MCcount, pip_MCcount, pim_MCcount, Kp_MCcount, Km_MCcount, g_MCcount;
int e_ind[BUFFER], p_ind[BUFFER], n_ind[BUFFER], pip_ind[BUFFER], pim_ind[BUFFER], Kp_ind[BUFFER], Km_ind[BUFFER], g_ind[BUFFER];

TLorentzVector p4_ele[BUFFER];
TLorentzVector p4_ele_raw[BUFFER];
TLorentzVector p4_prot[BUFFER];
TLorentzVector p4_neutr[BUFFER];
TLorentzVector p4_pip[BUFFER];
TLorentzVector p4_pim[BUFFER];
TLorentzVector p4_Kp[BUFFER];
TLorentzVector p4_Km[BUFFER];
TLorentzVector p4_phot[BUFFER];

float e_vx[BUFFER], e_vy[BUFFER], e_vz[BUFFER], e_beta[BUFFER], e_FTOF_sec[BUFFER], e_PCAL_sec[BUFFER];
float p_vx[BUFFER], p_vy[BUFFER], p_vz[BUFFER], p_beta[BUFFER], p_FTOF_sec[BUFFER], p_PCAL_sec[BUFFER];
float n_vx[BUFFER], n_vy[BUFFER], n_vz[BUFFER], n_beta[BUFFER];
float pip_vx[BUFFER], pip_vy[BUFFER], pip_vz[BUFFER], pip_beta[BUFFER], pip_FTOF_sec[BUFFER];
float pim_vx[BUFFER], pim_vy[BUFFER], pim_vz[BUFFER], pim_beta[BUFFER], pim_FTOF_sec[BUFFER];
float Kp_vx[BUFFER], Kp_vy[BUFFER], Kp_vz[BUFFER], Kp_beta[BUFFER], Kp_FTOF_sec[BUFFER];
float Km_vx[BUFFER], Km_vy[BUFFER], Km_vz[BUFFER], Km_beta[BUFFER], Km_FTOF_sec[BUFFER];
float g_vx[BUFFER], g_vy[BUFFER], g_vz[BUFFER], g_sec[BUFFER];

// cut variables

bool electron_sector_cut[BUFFER];


// vectors with components of the Lorentzvector to fill into the output tree
// --> vectors are used to optimize the file size

int eventnumber;
double beam_charge;
int helicity;
vector<double> p4_ele_px;
vector<double> p4_ele_py;
vector<double> p4_ele_pz;
vector<double> p4_ele_E;
vector<double> p4_prot_px;
vector<double> p4_prot_py;
vector<double> p4_prot_pz;
vector<double> p4_prot_E;
vector<double> p4_neutr_px;
vector<double> p4_neutr_py;
vector<double> p4_neutr_pz;
vector<double> p4_neutr_E;
vector<double> p4_pip_px;
vector<double> p4_pip_py;
vector<double> p4_pip_pz;
vector<double> p4_pip_E;
vector<double> p4_pim_px;
vector<double> p4_pim_py;
vector<double> p4_pim_pz;
vector<double> p4_pim_E;
vector<double> p4_Kp_px;
vector<double> p4_Kp_py;
vector<double> p4_Kp_pz;
vector<double> p4_Kp_E;
vector<double> p4_Km_px;
vector<double> p4_Km_py;
vector<double> p4_Km_pz;
vector<double> p4_Km_E;
vector<double> p4_phot_px;
vector<double> p4_phot_py;
vector<double> p4_phot_pz;
vector<double> p4_phot_E;
vector<double> p4_ele_chi2pid;
vector<double> p4_prot_chi2pid;
vector<double> p4_neutr_chi2pid;
vector<double> p4_pip_chi2pid;
vector<double> p4_pim_chi2pid;
vector<double> p4_Kp_chi2pid;
vector<double> p4_Km_chi2pid;
vector<double> p4_ele_dcx1;
vector<double> p4_prot_dcx1;
vector<double> p4_pip_dcx1; 
vector<double> p4_pim_dcx1; 
vector<double> p4_ele_dcy1;
vector<double> p4_prot_dcy1;
vector<double> p4_pip_dcy1; 
vector<double> p4_pim_dcy1; 
vector<double> p4_ele_dcx2;
vector<double> p4_pip_dcx2; 
vector<double> p4_pim_dcx2; 
vector<double> p4_ele_dcy2;
vector<double> p4_pip_dcy2; 
vector<double> p4_pim_dcy2; 
vector<double> p4_ele_dcx3;
vector<double> p4_pip_dcx3; 
vector<double> p4_pim_dcx3;
vector<double> p4_ele_dcy3;
vector<double> p4_pip_dcy3; 
vector<double> p4_pim_dcy3;
vector<double> p4_prot_dcz1;
vector<double> p4_ele_pcalv;
vector<double> p4_phot_pcalv;
vector<double> p4_ele_pcalw;
vector<double> p4_phot_pcalw;
vector<double> p4_ele_dcedge1;
vector<double> p4_ele_dcedge2;
vector<double> p4_ele_dcedge3;
vector<double> p4_prot_dcedge1;
vector<double> p4_prot_dcedge2;
vector<double> p4_prot_dcedge3;
vector<double> p4_pip_dcedge1;
vector<double> p4_pip_dcedge2;
vector<double> p4_pip_dcedge3;
vector<double> p4_pim_dcedge1;
vector<double> p4_pim_dcedge2;
vector<double> p4_pim_dcedge3;
vector<int> ele_det;
vector<int> prot_det;
vector<int> neutr_det;
vector<int> pip_det;
vector<int> pim_det;
vector<int> Kp_det;
vector<int> Km_det;
vector<int> phot_det;
vector<int> ele_sec;
vector<int> prot_sec;
vector<int> neutr_sec;
vector<int> pip_sec;
vector<int> pim_sec;
vector<int> Kp_sec;
vector<int> Km_sec;
vector<int> pip_RICHbest;
vector<int> pim_RICHbest;
vector<int> Kp_RICHbest;
vector<int> Km_RICHbest;
vector<int> phot_sec;

vector<double> p4_gen_ele_px;
vector<double> p4_gen_ele_py;
vector<double> p4_gen_ele_pz;
vector<double> p4_gen_ele_E;
vector<double> p4_gen_prot_px;
vector<double> p4_gen_prot_py;
vector<double> p4_gen_prot_pz;
vector<double> p4_gen_prot_E;
vector<double> p4_gen_neutr_px;
vector<double> p4_gen_neutr_py;
vector<double> p4_gen_neutr_pz;
vector<double> p4_gen_neutr_E;
vector<double> p4_gen_pip_px;
vector<double> p4_gen_pip_py;
vector<double> p4_gen_pip_pz;
vector<double> p4_gen_pip_E;
vector<double> p4_gen_pim_px;
vector<double> p4_gen_pim_py;
vector<double> p4_gen_pim_pz;
vector<double> p4_gen_pim_E;
vector<double> p4_gen_Kp_px;
vector<double> p4_gen_Kp_py;
vector<double> p4_gen_Kp_pz;
vector<double> p4_gen_Kp_E;
vector<double> p4_gen_Km_px;
vector<double> p4_gen_Km_py;
vector<double> p4_gen_Km_pz;
vector<double> p4_gen_Km_E;
vector<double> p4_gen_phot_px;
vector<double> p4_gen_phot_py;
vector<double> p4_gen_phot_pz;
vector<double> p4_gen_phot_E;
vector<double> p4_gen_pi0_px;
vector<double> p4_gen_pi0_py;
vector<double> p4_gen_pi0_pz;
vector<double> p4_gen_pi0_E;

vector<int> p4_gen_ele_motherpid;
vector<int> p4_gen_prot_motherpid;
vector<int> p4_gen_neutr_motherpid;
vector<int> p4_gen_pip_motherpid;
vector<int> p4_gen_pim_motherpid;
vector<int> p4_gen_Kp_motherpid;
vector<int> p4_gen_Km_motherpid;
vector<int> p4_gen_phot_motherpid;
vector<int> p4_gen_pi0_motherpid;

vector<double> p4_pip_ml;
vector<double> p4_pip_beta;
vector<double> p4_pip_LTCC_nphe;
vector<double> p4_pim_ml;
vector<double> p4_pim_beta;
vector<double> p4_pim_LTCC_nphe;
vector<double> p4_Kp_ml;
vector<double> p4_Kp_beta;
vector<double> p4_Kp_LTCC_nphe;
vector<double> p4_Km_ml;
vector<double> p4_Km_beta;
vector<double> p4_Km_LTCC_nphe;


// cut statistics:

long neg_part_count;
long pos_part_count;
long neut_part_count;

long Track_Quality_pass;

long FD_eid_default_PID_pass;
long FD_eid_charge_pass;
long FD_eid_CC_nphe_pass;
long FD_eid_EC_outer_vs_EC_inner_pass;
long FD_eid_EC_sampling_fraction_pass;
long FD_eid_EC_hit_position_fiducial_pass;
long FD_eid_DC_hit_position_region1_fiducial_pass;
long FD_eid_DC_hit_position_region2_fiducial_pass;
long FD_eid_DC_hit_position_region3_fiducial_pass;
long FD_eid_DC_z_vertex_pass;
long FD_eid_all_pass;

long FD_protid_default_PID_pass;
long FD_protid_charge_pass;
long FD_protid_DC_hit_position_region1_fiducial_pass;
long FD_protid_DC_hit_position_region2_fiducial_pass;
long FD_protid_DC_hit_position_region3_fiducial_pass;
long FD_protid_beta_pass;
long FD_protid_delta_beta_pass;
long FD_protid_tofmass_pass;
long FD_protid_maximum_probability_pass;
long FD_protid_delta_vz_pass;
long FD_protid_all_pass;

long FD_neutrid_default_PID_pass;
long FD_neutrid_charge_pass;
long FD_neutrid_beta_pass;
long FD_neutrid_delta_beta_pass;
long FD_neutrid_tofmass_pass;
long FD_neutrid_delta_vz_pass;
long FD_neutrid_all_pass;

long FD_pipid_default_PID_pass;
long FD_pipid_charge_pass;
long FD_pipid_DC_hit_position_region1_fiducial_pass;
long FD_pipid_DC_hit_position_region2_fiducial_pass;
long FD_pipid_DC_hit_position_region3_fiducial_pass;
long FD_pipid_beta_pass;
long FD_pipid_delta_beta_pass;
long FD_pipid_tofmass_pass;
long FD_pipid_maximum_probability_pass;
long FD_pipid_delta_vz_pass;
long FD_pipid_all_pass;

long FD_pimid_default_PID_pass;
long FD_pimid_charge_pass;
long FD_pimid_DC_hit_position_region1_fiducial_pass;
long FD_pimid_DC_hit_position_region2_fiducial_pass;
long FD_pimid_DC_hit_position_region3_fiducial_pass;
long FD_pimid_beta_pass;
long FD_pimid_delta_beta_pass;
long FD_pimid_tofmass_pass;
long FD_pimid_maximum_probability_pass;
long FD_pimid_delta_vz_pass;
long FD_pimid_all_pass;

long FD_Kpid_default_PID_pass;
long FD_Kpid_charge_pass;
long FD_Kpid_DC_hit_position_region1_fiducial_pass;
long FD_Kpid_DC_hit_position_region2_fiducial_pass;
long FD_Kpid_DC_hit_position_region3_fiducial_pass;
long FD_Kpid_beta_pass;
long FD_Kpid_delta_beta_pass;
long FD_Kpid_tofmass_pass;
long FD_Kpid_maximum_probability_pass;
long FD_Kpid_delta_vz_pass;
long FD_Kpid_all_pass;

long FD_Kmid_default_PID_pass;
long FD_Kmid_charge_pass;
long FD_Kmid_DC_hit_position_region1_fiducial_pass;
long FD_Kmid_DC_hit_position_region2_fiducial_pass;
long FD_Kmid_DC_hit_position_region3_fiducial_pass;
long FD_Kmid_beta_pass;
long FD_Kmid_delta_beta_pass;
long FD_Kmid_tofmass_pass;
long FD_Kmid_maximum_probability_pass;
long FD_Kmid_delta_vz_pass;
long FD_Kmid_all_pass;

long FD_photid_default_PID_pass;
long FD_photid_charge_pass;
long FD_photid_beta_pass;
long FD_photid_EC_sampling_fraction_pass;
long FD_photid_EC_hit_position_fiducial_pass;
long FD_photid_all_pass;

/// FT

long FT_eid_charge_pass;
long FT_eid_PID_pass;
long FT_eid_FTCAL_fiducial_pass;
long FT_eid_FTTRK_fiducial_pass;
long FT_eid_FTHODO_fiducial_pass;
long FT_eid_energy_vs_radius_pass;
long FT_eid_all_pass;

long FT_photid_charge_pass;
long FT_photid_PID_pass;
long FT_photid_FTCAL_fiducial_pass;
long FT_photid_beta_pass;
long FT_photid_all_pass;

/// CD

long CD_protid_default_PID_pass;
long CD_protid_charge_pass;
long CD_protid_beta_pass;
long CD_protid_maximum_probability_pass;
long CD_protid_delta_vz_pass;
long CD_protid_all_pass;

long CD_neutrid_default_PID_pass;
long CD_neutrid_charge_pass;
long CD_neutrid_beta_pass;
long CD_neutrid_maximum_probability_pass;
long CD_neutrid_delta_vz_pass;
long CD_neutrid_all_pass;

long CD_pipid_default_PID_pass;
long CD_pipid_charge_pass;
long CD_pipid_beta_pass;
long CD_pipid_maximum_probability_pass;
long CD_pipid_delta_vz_pass;
long CD_pipid_all_pass;

long CD_pimid_default_PID_pass;
long CD_pimid_charge_pass;
long CD_pimid_beta_pass;
long CD_pimid_maximum_probability_pass;
long CD_pimid_delta_vz_pass;
long CD_pimid_all_pass;

long CD_Kpid_default_PID_pass;
long CD_Kpid_charge_pass;
long CD_Kpid_beta_pass;
long CD_Kpid_maximum_probability_pass;
long CD_Kpid_delta_vz_pass;
long CD_Kpid_all_pass;

long CD_Kmid_default_PID_pass;
long CD_Kmid_charge_pass;
long CD_Kmid_beta_pass;
long CD_Kmid_maximum_probability_pass;
long CD_Kmid_delta_vz_pass;
long CD_Kmid_all_pass;

// cut selectors variables

bool Track_Quality_check[BUFFER];

/// /////////////////////////////////////////////////////////////
/// Forward detector

bool FD_eid_default_PID_check[BUFFER];
bool FD_eid_charge_check[BUFFER];
bool FD_eid_CC_nphe_check[BUFFER];
bool FD_eid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_eid_EC_sampling_fraction_check[BUFFER];
bool FD_eid_EC_hit_position_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_eid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_eid_DC_z_vertex_check[BUFFER];
bool FD_eid_all_check[BUFFER];

bool FD_protid_default_PID_check[BUFFER];
bool FD_protid_charge_check[BUFFER];
bool FD_protid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_protid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_protid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_protid_beta_check[BUFFER];
bool FD_protid_delta_beta_check[BUFFER];
bool FD_protid_tofmass_check[BUFFER];
bool FD_protid_maximum_probability_check[BUFFER];
bool FD_protid_delta_vz_check[BUFFER];
bool FD_protid_all_check[BUFFER];

bool FD_neutrid_default_PID_check[BUFFER];
bool FD_neutrid_charge_check[BUFFER];
bool FD_neutrid_beta_check[BUFFER];
bool FD_neutrid_delta_beta_check[BUFFER];
bool FD_neutrid_tofmass_check[BUFFER];
bool FD_neutrid_delta_vz_check[BUFFER];
bool FD_neutrid_all_check[BUFFER];

bool FD_pipid_default_PID_check[BUFFER];
bool FD_pipid_charge_check[BUFFER];
bool FD_pipid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_pipid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_pipid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_pipid_beta_check[BUFFER];
bool FD_pipid_delta_beta_check[BUFFER];
bool FD_pipid_tofmass_check[BUFFER];
bool FD_pipid_maximum_probability_check[BUFFER];
bool FD_pipid_delta_vz_check[BUFFER];
bool FD_pipid_all_check[BUFFER];

bool FD_pimid_default_PID_check[BUFFER];
bool FD_pimid_charge_check[BUFFER];
bool FD_pimid_ele_reject_check[BUFFER];
bool FD_pimid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_pimid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_pimid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_pimid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_pimid_beta_check[BUFFER];
bool FD_pimid_delta_beta_check[BUFFER];
bool FD_pimid_tofmass_check[BUFFER];
bool FD_pimid_maximum_probability_check[BUFFER];
bool FD_pimid_delta_vz_check[BUFFER];
bool FD_pimid_all_check[BUFFER];

bool FD_Kpid_default_PID_check[BUFFER];
bool FD_Kpid_charge_check[BUFFER];
bool FD_Kpid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_Kpid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_Kpid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_Kpid_beta_check[BUFFER];
bool FD_Kpid_delta_beta_check[BUFFER];
bool FD_Kpid_tofmass_check[BUFFER];
bool FD_Kpid_maximum_probability_check[BUFFER];
bool FD_Kpid_delta_vz_check[BUFFER];
bool FD_Kpid_all_check[BUFFER];

bool FD_Kmid_default_PID_check[BUFFER];
bool FD_Kmid_charge_check[BUFFER];
bool FD_Kmid_ele_reject_check[BUFFER];
bool FD_Kmid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_Kmid_DC_hit_position_region1_fiducial_check[BUFFER];
bool FD_Kmid_DC_hit_position_region2_fiducial_check[BUFFER];
bool FD_Kmid_DC_hit_position_region3_fiducial_check[BUFFER];
bool FD_Kmid_beta_check[BUFFER];
bool FD_Kmid_delta_beta_check[BUFFER];
bool FD_Kmid_tofmass_check[BUFFER];
bool FD_Kmid_maximum_probability_check[BUFFER];
bool FD_Kmid_delta_vz_check[BUFFER];
bool FD_Kmid_all_check[BUFFER];

bool FD_photid_default_PID_check[BUFFER];
bool FD_photid_charge_check[BUFFER];
bool FD_photid_beta_check[BUFFER];
bool FD_photid_EC_sampling_fraction_check[BUFFER];
bool FD_photid_EC_outer_vs_EC_inner_check[BUFFER];
bool FD_photid_EC_hit_position_fiducial_check[BUFFER];
bool FD_photid_all_check[BUFFER];

/// //////////////////////////////////////////////////////////
/// FT

bool FT_eid_charge_check[BUFFER];
bool FT_eid_PID_check[BUFFER];
bool FT_eid_FTCAL_fiducial_check[BUFFER];
bool FT_eid_FTTRK_fiducial_check[BUFFER];
bool FT_eid_FTHODO_fiducial_check[BUFFER];
bool FT_eid_energy_vs_radius_check[BUFFER];
bool FT_eid_all_check[BUFFER];

bool FT_photid_charge_check[BUFFER];
bool FT_photid_PID_check[BUFFER];
bool FT_photid_FTCAL_fiducial_check[BUFFER];
bool FT_photid_beta_check[BUFFER];
bool FT_photid_all_check[BUFFER];

/// //////////////////////////////////////////////////////////////
/// Central detector

bool CD_protid_default_PID_check[BUFFER];
bool CD_protid_charge_check[BUFFER];
bool CD_protid_beta_check[BUFFER];
bool CD_protid_maximum_probability_check[BUFFER];
bool CD_protid_delta_vz_check[BUFFER];
bool CD_protid_all_check[BUFFER];

bool CD_neutrid_default_PID_check[BUFFER];
bool CD_neutrid_charge_check[BUFFER];
bool CD_neutrid_beta_check[BUFFER];
bool CD_neutrid_maximum_probability_check[BUFFER];
bool CD_neutrid_delta_vz_check[BUFFER];
bool CD_neutrid_all_check[BUFFER];

bool CD_pipid_default_PID_check[BUFFER];
bool CD_pipid_charge_check[BUFFER];
bool CD_pipid_beta_check[BUFFER];
bool CD_pipid_maximum_probability_check[BUFFER];
bool CD_pipid_delta_vz_check[BUFFER];
bool CD_pipid_all_check[BUFFER];

bool CD_pimid_default_PID_check[BUFFER];
bool CD_pimid_charge_check[BUFFER];
bool CD_pimid_beta_check[BUFFER];
bool CD_pimid_maximum_probability_check[BUFFER];
bool CD_pimid_delta_vz_check[BUFFER];
bool CD_pimid_all_check[BUFFER];

bool CD_Kpid_default_PID_check[BUFFER];
bool CD_Kpid_charge_check[BUFFER];
bool CD_Kpid_beta_check[BUFFER];
bool CD_Kpid_maximum_probability_check[BUFFER];
bool CD_Kpid_delta_vz_check[BUFFER];
bool CD_Kpid_all_check[BUFFER];

bool CD_Kmid_default_PID_check[BUFFER];
bool CD_Kmid_charge_check[BUFFER];
bool CD_Kmid_beta_check[BUFFER];
bool CD_Kmid_maximum_probability_check[BUFFER];
bool CD_Kmid_delta_vz_check[BUFFER];
bool CD_Kmid_all_check[BUFFER];

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
/// PID histograms:

/// a) electron ID:

const int FD_eid_cuts = 16;

// CC cuts

TH2F *hist_HTCC_theta_vs_phi[FD_eid_cuts];
TH1F *hist_HTCC_nphe[FD_eid_cuts];
TH2F *hist_HTCC_nphe_vs_sampling_fraction[FD_eid_cuts];

// EC cuts

TH2F *hist_EC_PCAL_vs_EC_ECAL[FD_eid_cuts];
TH2F *hist_EC_outer_vs_EC_inner[FD_eid_cuts];

TH2F *hist_EC_total_sampling_fraction_sec1[FD_eid_cuts];
TH2F *hist_EC_total_sampling_fraction_sec2[FD_eid_cuts];
TH2F *hist_EC_total_sampling_fraction_sec3[FD_eid_cuts];
TH2F *hist_EC_total_sampling_fraction_sec4[FD_eid_cuts];
TH2F *hist_EC_total_sampling_fraction_sec5[FD_eid_cuts];
TH2F *hist_EC_total_sampling_fraction_sec6[FD_eid_cuts];

TH2F *hist_EC_PCAL_sampling_fraction_sec1[FD_eid_cuts];
TH2F *hist_EC_PCAL_sampling_fraction_sec2[FD_eid_cuts];
TH2F *hist_EC_PCAL_sampling_fraction_sec3[FD_eid_cuts];
TH2F *hist_EC_PCAL_sampling_fraction_sec4[FD_eid_cuts];
TH2F *hist_EC_PCAL_sampling_fraction_sec5[FD_eid_cuts];
TH2F *hist_EC_PCAL_sampling_fraction_sec6[FD_eid_cuts];

TH2F *hist_EC_ECAL_sampling_fraction_sec1[FD_eid_cuts];
TH2F *hist_EC_ECAL_sampling_fraction_sec2[FD_eid_cuts];
TH2F *hist_EC_ECAL_sampling_fraction_sec3[FD_eid_cuts];
TH2F *hist_EC_ECAL_sampling_fraction_sec4[FD_eid_cuts];
TH2F *hist_EC_ECAL_sampling_fraction_sec5[FD_eid_cuts];
TH2F *hist_EC_ECAL_sampling_fraction_sec6[FD_eid_cuts];

TH2F *hist_EC_PCAL_hit_position[FD_eid_cuts];
TH2F *hist_EC_inner_hit_position[FD_eid_cuts];
TH2F *hist_EC_outer_hit_position[FD_eid_cuts];

// DC cuts

TH2F *hist_DC_hit_position_region1[FD_eid_cuts];
TH2F *hist_DC_hit_position_region2[FD_eid_cuts];
TH2F *hist_DC_hit_position_region3[FD_eid_cuts];

TH1F *hist_DC_z_vertex_sec1[FD_eid_cuts];
TH1F *hist_DC_z_vertex_sec2[FD_eid_cuts];
TH1F *hist_DC_z_vertex_sec3[FD_eid_cuts];
TH1F *hist_DC_z_vertex_sec4[FD_eid_cuts];
TH1F *hist_DC_z_vertex_sec5[FD_eid_cuts];
TH1F *hist_DC_z_vertex_sec6[FD_eid_cuts];

/// b) FTOF + additional cuts for charged hadrons:

const int FD_hid_count = 60;

TH2F *hist_DC_hit_position_region1_hadron[FD_hid_count];
TH2F *hist_DC_hit_position_region2_hadron[FD_hid_count];
TH2F *hist_DC_hit_position_region3_hadron[FD_hid_count];
TH2F *hist_EC_outer_vs_EC_inner_hadron[FD_hid_count];

TH2F *hist_beta_vs_p[FD_hid_count];
TH2F *hist_beta_vs_p_sec1[FD_hid_count];
TH2F *hist_beta_vs_p_sec2[FD_hid_count];
TH2F *hist_beta_vs_p_sec3[FD_hid_count];
TH2F *hist_beta_vs_p_sec4[FD_hid_count];
TH2F *hist_beta_vs_p_sec5[FD_hid_count];
TH2F *hist_beta_vs_p_sec6[FD_hid_count];

TH2F *hist_delta_beta_vs_p[FD_hid_count];
TH1F *hist_delta_beta[FD_hid_count];
TH2F *hist_tofmass_vs_p[FD_hid_count];
TH1F *hist_tofmass[FD_hid_count];
TH1F *hist_delta_vz[FD_hid_count];

/// CTOF + additional cuts for charged hadrons:

const int FD_hid_CD_count = 60;

TH2F *hist_CD_beta_vs_p[FD_hid_CD_count];
TH1F *hist_CD_delta_vz[FD_hid_CD_count];

/// c) photon ID:

const int FD_photid_count = 8;

TH2F *hist_beta_vs_p_phot[FD_photid_count];
TH1F *hist_beta_phot[FD_photid_count];
TH2F *hist_EC_sampling_fraction_phot[FD_photid_count];
TH2F *hist_EC_PCAL_vs_EC_ECAL_phot[FD_photid_count];
TH2F *hist_EC_PCAL_hit_position_phot[FD_photid_count];

////////////////////////////////////////////////////////////////////
///  FT plots

const int FT_pid_count = 20;

TH2F *hist_FT_FTCAL_energy_vs_radius[FT_pid_count];
TH2F *hist_FT_FTCAL_hit_position[FT_pid_count];
TH2F *hist_FT_FTTRK_hit_position[FT_pid_count];
TH2F *hist_FT_FTHODO_hit_position[FT_pid_count];
TH1F *hist_FT_beta[FT_pid_count];

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// define functions

/// electron ID:

// CC cuts

TH2F *create_hist_HTCC_theta_vs_phi(int cutnum);
TH1F *create_hist_HTCC_nphe(int cutnum);
TH2F *create_hist_HTCC_nphe_vs_sampling_fraction(int cutnum);

// EC cuts

TH2F *create_hist_EC_PCAL_vs_EC_ECAL(int cutnum);
TH2F *create_hist_EC_outer_vs_EC_inner(int cutnum);

TH2F *create_hist_EC_total_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_total_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_PCAL_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_PCAL_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_ECAL_sampling_fraction_sec1(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec2(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec3(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec4(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec5(int cutnum);
TH2F *create_hist_EC_ECAL_sampling_fraction_sec6(int cutnum);

TH2F *create_hist_EC_PCAL_hit_position(int cutnum);
TH2F *create_hist_EC_inner_hit_position(int cutnum);
TH2F *create_hist_EC_outer_hit_position(int cutnum);

// DC cuts

TH2F *create_hist_DC_hit_position_region1(int cutnum);
TH2F *create_hist_DC_hit_position_region2(int cutnum);
TH2F *create_hist_DC_hit_position_region3(int cutnum);

TH1F *create_hist_DC_z_vertex_sec1(int cutnum);
TH1F *create_hist_DC_z_vertex_sec2(int cutnum);
TH1F *create_hist_DC_z_vertex_sec3(int cutnum);
TH1F *create_hist_DC_z_vertex_sec4(int cutnum);
TH1F *create_hist_DC_z_vertex_sec5(int cutnum);
TH1F *create_hist_DC_z_vertex_sec6(int cutnum);

// TOF + others for charged hadrons

TH2F *create_DC_hit_position_region1_hadron(int cutnum);
TH2F *create_DC_hit_position_region2_hadron(int cutnum);
TH2F *create_DC_hit_position_region3_hadron(int cutnum);
TH2F *create_hist_EC_outer_vs_EC_inner_hadron(int cutnum);

TH2F *create_hist_beta_vs_p(int cutnum);
TH2F *create_hist_beta_vs_p_sec1(int cutnum);
TH2F *create_hist_beta_vs_p_sec2(int cutnum);
TH2F *create_hist_beta_vs_p_sec3(int cutnum);
TH2F *create_hist_beta_vs_p_sec4(int cutnum);
TH2F *create_hist_beta_vs_p_sec5(int cutnum);
TH2F *create_hist_beta_vs_p_sec6(int cutnum);

TH2F *create_hist_delta_beta_vs_p(int cutnum);
TH1F *create_hist_delta_beta(int cutnum);
TH2F *create_hist_tofmass_vs_p(int cutnum);
TH1F *create_hist_tofmass(int cutnum);
TH1F *create_hist_delta_vz(int cutnum);

// CD charged hadrons

TH2F *create_hist_CD_beta_vs_p(int cutnum);
TH1F *create_hist_CD_delta_vz(int cutnum);

// TOF for photons

TH2F *create_hist_beta_vs_p_phot(int cutnum);
TH1F *create_hist_beta_phot(int cutnum);
TH2F *create_hist_EC_sampling_fraction_phot(int cutnum);
TH2F *create_hist_EC_PCAL_vs_EC_ECAL_phot(int cutnum);
TH2F *create_hist_EC_PCAL_hit_position_phot(int cutnum);

////////////////////////////////////////////////////////////////
// FT plots

TH2F *create_hist_FT_FTCAL_energy_vs_radius(int cutnum);
TH2F *create_hist_FT_FTCAL_hit_position(int cutnum);
TH2F *create_hist_FT_FTTRK_hit_position(int cutnum);
TH2F *create_hist_FT_FTHODO_hit_position(int cutnum);
TH1F *create_hist_FT_beta(int cutnum);

// get event properties from the EVNT bank

void get_event_properties(clas12::event_ptr event, clas12::mcevt_ptr mcevent);

// raw particle assignment

void assign_particles(std::vector<region_part_ptr>& particles, mcpar_ptr mcparticles);

// particle selection:

void select_electron(int run);
void select_proton(int run);
void select_neutron(int run);
void select_pip(int run);
void select_pim(int run);
void select_Kplus(int run);
void select_Kminus(int run);
void select_photon(int run);

void write_gen(int run);

// momentum correction:

TLorentzVector correct_lepton_negative(double thetaeld, double phield, double pel, int secte);
TVector3 correct_hadron_positive(double thetahd, double phihd, double ph, int secth);
TVector3 correct_hadron_negative(double thetahd, double phihd, double ph, int secth);

void correct_electron(void);
void correct_proton(void);
void correct_neutron(void);
void correct_pip(void);
void correct_pim(void);
void correct_Kplus(void);
void correct_Kminus(void);
void correct_photon(void);

// fill output variables:

void fill_output_vector_electron(void);
void fill_output_vector_proton(void);
void fill_output_vector_neutron(void);
void fill_output_vector_pip(void);
void fill_output_vector_pim(void);
void fill_output_vector_Kplus(void);
void fill_output_vector_Kminus(void);
void fill_output_vector_photon(void);
void fill_output_vector_MC(void);

// define cuts:

bool basic_FTOF_cut(int j);

// PID checks:

bool Track_Quality_cut(int j);

bool ele_default_PID_cut(int j);
bool ele_charge_cut(int j);
bool CC_nphe_cut(int j);
bool EC_outer_vs_EC_inner_cut(int j);
bool EC_sampling_fraction_cut(int j);
bool EC_hit_position_fiducial_cut(int j);
bool EC_hit_position_fiducial_cut_homogeneous(int j);
bool DC_fiducial_cut_XY(int j, int region);
bool DC_fiducial_cut_theta_phi(int j, int region);
bool DC_fiducial_cut_edge(int j, int region);
bool DC_z_vertex_cut(int j);

bool prot_default_PID_cut(int j);
bool prot_charge_cut(int j);
bool prot_delta_vz_cut(int j);

bool neutr_default_PID_cut(int j);
bool neutr_charge_cut(int j);
bool neutr_beta_cut(int j, int run);
bool neutr_delta_beta_cut(int j, int run);
bool neutr_tofmass_cut(int j, int run);
bool neutr_delta_vz_cut(int j);

bool pip_default_PID_cut(int j);
bool pip_charge_cut(int j);
bool pip_delta_vz_cut(int j);

bool pim_default_PID_cut(int j);
bool pim_charge_cut(int j);
bool pim_ele_reject_cut(int j);
bool pim_EC_outer_vs_EC_inner_cut(int j);
bool pim_delta_vz_cut(int j);

bool Kp_default_PID_cut(int j);
bool Kp_charge_cut(int j);
bool Kp_delta_vz_cut(int j);

bool Km_default_PID_cut(int j);
bool Km_charge_cut(int j);
bool Km_ele_reject_cut(int j);
bool Km_EC_outer_vs_EC_inner_cut(int j);
bool Km_delta_vz_cut(int j);

bool maximum_probability_cut(int j, int hypothesis, double conflvl, double anticonflvl, int run);

bool phot_default_PID_cut(int j);
bool phot_charge_cut(int j);
bool phot_beta_cut(int j, int run);
bool phot_EC_sampling_fraction_cut(int j);
bool phot_EC_outer_vs_EC_inner_cut(int j);
bool phot_EC_hit_position_fiducial_cut(int j);

// FT

bool FT_eid_charge_cut(int j);
bool FT_eid_PID_cut(int j);
bool FT_eid_FTCAL_fiducial_cut(int j);
bool FT_eid_FTTRK_fiducial_cut(int j);
bool FT_eid_FTHODO_fiducial_cut(int j);
bool FT_eid_energy_vs_radius_cut(int j);

bool FT_photid_charge_cut(int j);
bool FT_photid_PID_cut(int j);
bool FT_photid_FTCAL_fiducial_cut(int j);
bool FT_photid_beta_cut(int j, int run);

// CD

bool CD_prot_default_PID_cut(int j);
bool CD_prot_charge_cut(int j);
bool CD_prot_beta_cut(int j, int run);
bool CD_prot_delta_vz_cut(int j);

bool CD_neutr_default_PID_cut(int j);
bool CD_neutr_charge_cut(int j);
bool CD_neutr_beta_cut(int j, int run);
bool CD_neutr_delta_vz_cut(int j);

bool CD_pip_default_PID_cut(int j);
bool CD_pip_charge_cut(int j);
bool CD_pip_beta_cut(int j, int run);
bool CD_pip_delta_vz_cut(int j);

bool CD_pim_default_PID_cut(int j);
bool CD_pim_charge_cut(int j);
bool CD_pim_beta_cut(int j, int run);
bool CD_pim_delta_vz_cut(int j);

bool CD_Kp_default_PID_cut(int j);
bool CD_Kp_charge_cut(int j);
bool CD_Kp_beta_cut(int j, int run);
bool CD_Kp_delta_vz_cut(int j);

bool CD_Km_default_PID_cut(int j);
bool CD_Km_charge_cut(int j);
bool CD_Km_beta_cut(int j, int run);
bool CD_Km_delta_vz_cut(int j);

int determineSector(int i);

bool good_sc_paddle(int j);

// additional calculation functions:

double GetTheta(int j);
double GetPhi(int j);
TVector3 GetUVWVector(int j);
double Getdvz(int j);
double Beta_charged(int j, int run);
double Beta_charged_central(int j, int run);
double Beta_charged_FT(int j, int run);
double Beta_neutral(int j, int run);
double Beta_neutral_FT(int j, int run);
double Get_Starttime(int j, int run);


double ml_value(int i);

/// /////////////////////////////////////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// main
///

Int_t filter_clas12_hipo(Char_t *inFile, Char_t *outputfile, int run)
{

    clas12::clas12reader c12(inFile); // open file

    clas12databases db;
    c12.connectDataBases(&db);

    cout << "Initalize the input file ... " << endl;

    Char_t tmpstr[80];
    Double_t fraction;

    /// /////////////////////////////////////////////////////////////////////////////
    ///  create output-txtfiles:
    /// /////////////////////////////////////////////////////////////////////////////

    // ofstream outputFile_electrons_elastic("output/electrons_elastic.txt");
    // ofstream outputFile_ep_elastic("output/ep_elastic.txt");
    // ofstream outputFile_epip_missing_neutron("output/e_pip_X.txt");

    /// /////////////////////////////////////////////////////////////////////////////
    ///  create output-file and tree for saving histograms:
    /// /////////////////////////////////////////////////////////////////////////////

    out = new TFile(outputfile, "RECREATE");

    TTree out_tree("out_tree", "out_tree");
    out_tree.Branch("eventnumber", &eventnumber);
    out_tree.Branch("helicity", &helicity);
    out_tree.Branch("beam_charge", &beam_charge);
    out_tree.Branch("p4_ele_px", &p4_ele_px);
    out_tree.Branch("p4_ele_py", &p4_ele_py);
    out_tree.Branch("p4_ele_pz", &p4_ele_pz);
    out_tree.Branch("p4_ele_E", &p4_ele_E);
    out_tree.Branch("p4_prot_px", &p4_prot_px);
    out_tree.Branch("p4_prot_py", &p4_prot_py);
    out_tree.Branch("p4_prot_pz", &p4_prot_pz);
    out_tree.Branch("p4_prot_E", &p4_prot_E);
    out_tree.Branch("p4_neutr_px", &p4_neutr_px);
    out_tree.Branch("p4_neutr_py", &p4_neutr_py);
    out_tree.Branch("p4_neutr_pz", &p4_neutr_pz);
    out_tree.Branch("p4_neutr_E", &p4_neutr_E);
    out_tree.Branch("p4_pip_px", &p4_pip_px);
    out_tree.Branch("p4_pip_py", &p4_pip_py);
    out_tree.Branch("p4_pip_pz", &p4_pip_pz);
    out_tree.Branch("p4_pip_E", &p4_pip_E);
    out_tree.Branch("p4_pim_px", &p4_pim_px);
    out_tree.Branch("p4_pim_py", &p4_pim_py);
    out_tree.Branch("p4_pim_pz", &p4_pim_pz);
    out_tree.Branch("p4_pim_E", &p4_pim_E);
    out_tree.Branch("p4_Kp_px", &p4_Kp_px);
    out_tree.Branch("p4_Kp_py", &p4_Kp_py);
    out_tree.Branch("p4_Kp_pz", &p4_Kp_pz);
    out_tree.Branch("p4_Kp_E", &p4_Kp_E);
    out_tree.Branch("p4_Km_px", &p4_Km_px);
    out_tree.Branch("p4_Km_py", &p4_Km_py);
    out_tree.Branch("p4_Km_pz", &p4_Km_pz);
    out_tree.Branch("p4_Km_E", &p4_Km_E);
    out_tree.Branch("p4_phot_px", &p4_phot_px);
    out_tree.Branch("p4_phot_py", &p4_phot_py);
    out_tree.Branch("p4_phot_pz", &p4_phot_pz);
    out_tree.Branch("p4_phot_E", &p4_phot_E);
    out_tree.Branch("ele_det", &ele_det);
    out_tree.Branch("prot_det", &prot_det);
    out_tree.Branch("neutr_det", &neutr_det);
    out_tree.Branch("pip_det", &pip_det);
    out_tree.Branch("pim_det", &pim_det);
    out_tree.Branch("Kp_det", &Kp_det);
    out_tree.Branch("Km_det", &Km_det);
    out_tree.Branch("phot_det", &phot_det);
    out_tree.Branch("ele_sec", &ele_sec);
    out_tree.Branch("prot_sec", &prot_sec);
    out_tree.Branch("neutr_sec", &neutr_sec);
    out_tree.Branch("pip_sec", &pip_sec);
    out_tree.Branch("pim_sec", &pim_sec);
    out_tree.Branch("Kp_sec", &Kp_sec);
    out_tree.Branch("Km_sec", &Km_sec);
    out_tree.Branch("phot_sec", &phot_sec);
    out_tree.Branch("p4_ele_chi2pid", &p4_ele_chi2pid);
    out_tree.Branch("p4_prot_chi2pid", &p4_prot_chi2pid);
    out_tree.Branch("p4_neutr_chi2pid", &p4_neutr_chi2pid);
    out_tree.Branch("p4_pip_chi2pid", &p4_pip_chi2pid);
    out_tree.Branch("p4_pim_chi2pid", &p4_pim_chi2pid);
    out_tree.Branch("p4_Kp_chi2pid", &p4_Kp_chi2pid);
    out_tree.Branch("p4_Km_chi2pid", &p4_Km_chi2pid);
    if(writebeta == true){
      out_tree.Branch("p4_Kp_beta", &p4_Kp_beta);
      out_tree.Branch("p4_pip_beta", &p4_pip_beta);
      out_tree.Branch("p4_Km_beta", &p4_Km_beta);
      out_tree.Branch("p4_pim_beta", &p4_pim_beta);
    }
    if(writeML == true){
      out_tree.Branch("p4_pip_ml", &p4_pip_ml);
      out_tree.Branch("p4_pim_ml", &p4_pim_ml);
      out_tree.Branch("p4_Kp_ml", &p4_Kp_ml);
      out_tree.Branch("p4_Km_ml", &p4_Km_ml);
    }
    if(writeLTCC){    
      out_tree.Branch("p4_Kp_LTCC_nphe", &p4_Kp_LTCC_nphe);
      out_tree.Branch("p4_pip_LTCC_nphe", &p4_pip_LTCC_nphe);
      out_tree.Branch("p4_Km_LTCC_nphe", &p4_Km_LTCC_nphe);
      out_tree.Branch("p4_pim_LTCC_nphe", &p4_pim_LTCC_nphe);
    }
    if(userich == true){
      out_tree.Branch("pip_RICHbest", &pip_RICHbest);
      out_tree.Branch("pim_RICHbest", &pim_RICHbest);
      out_tree.Branch("Kp_RICHbest", &Kp_RICHbest);
      out_tree.Branch("Km_RICHbest", &Km_RICHbest);
    }
    if(fill_coord_ele == true){
      out_tree.Branch("p4_ele_dcx1", &p4_ele_dcx1);
      out_tree.Branch("p4_ele_dcy1", &p4_ele_dcy1);
      out_tree.Branch("p4_ele_dcx2", &p4_ele_dcx2);
      out_tree.Branch("p4_ele_dcy2", &p4_ele_dcy2);
      out_tree.Branch("p4_ele_dcx3", &p4_ele_dcx3);
      out_tree.Branch("p4_ele_dcy3", &p4_ele_dcy3);
      out_tree.Branch("p4_ele_pcalv", &p4_ele_pcalv);
      out_tree.Branch("p4_ele_pcalw", &p4_ele_pcalw);
    }
    if(fill_coord_prot == true){
      out_tree.Branch("p4_prot_dcx1", &p4_prot_dcx1);
      out_tree.Branch("p4_prot_dcy1", &p4_prot_dcy1);
      out_tree.Branch("p4_prot_dcz1", &p4_prot_dcz1);
    }      
    if(fill_coord_pip == true){          
      out_tree.Branch("p4_pip_dcx1", &p4_pip_dcx1);
      out_tree.Branch("p4_pip_dcy1", &p4_pip_dcy1);
      out_tree.Branch("p4_pip_dcx2", &p4_pip_dcx2);
      out_tree.Branch("p4_pip_dcy2", &p4_pip_dcy2);
      out_tree.Branch("p4_pip_dcx3", &p4_pip_dcx3);
      out_tree.Branch("p4_pip_dcy3", &p4_pip_dcy3);
    }
    if(fill_coord_pim == true){
      out_tree.Branch("p4_pim_dcx1", &p4_pim_dcx1);
      out_tree.Branch("p4_pim_dcy1", &p4_pim_dcy1);
      out_tree.Branch("p4_pim_dcx2", &p4_pim_dcx2);
      out_tree.Branch("p4_pim_dcy2", &p4_pim_dcy2);
      out_tree.Branch("p4_pim_dcx3", &p4_pim_dcx3);
      out_tree.Branch("p4_pim_dcy3", &p4_pim_dcy3);
    }
    if(fill_coord_phot == true){
      out_tree.Branch("p4_phot_pcalv", &p4_phot_pcalv); 
      out_tree.Branch("p4_phot_pcalw", &p4_phot_pcalw); 
    }
    if(fill_edge_ele == true){
      out_tree.Branch("p4_ele_dcedge1", &p4_ele_dcedge1);
      out_tree.Branch("p4_ele_dcedge2", &p4_ele_dcedge2);
      out_tree.Branch("p4_ele_dcedge3", &p4_ele_dcedge3);
    }
    if(fill_edge_prot == true){
      out_tree.Branch("p4_prot_dcedge1", &p4_prot_dcedge1);
      out_tree.Branch("p4_prot_dcedge2", &p4_prot_dcedge2);
      out_tree.Branch("p4_prot_dcedge3", &p4_prot_dcedge3);
    }      
    if(fill_edge_pip == true){          
      out_tree.Branch("p4_pip_dcedge1", &p4_pip_dcedge1);
      out_tree.Branch("p4_pip_dcedge2", &p4_pip_dcedge2);
      out_tree.Branch("p4_pip_dcedge3", &p4_pip_dcedge3);
    }
    if(fill_edge_pim == true){
      out_tree.Branch("p4_pim_dcedge1", &p4_pim_dcedge1);
      out_tree.Branch("p4_pim_dcedge2", &p4_pim_dcedge2);
      out_tree.Branch("p4_pim_dcedge3", &p4_pim_dcedge3);
    }
    if (simulation == true){
        out_tree.Branch("gen_event_helicity", &MC_helicity);
        out_tree.Branch("gen_event_npart", &MC_Npart);
        out_tree.Branch("gen_event_ebeam", &MC_Ebeam);
        out_tree.Branch("gen_event_weight", &MC_weight);
        out_tree.Branch("p4_gen_ele_px", &p4_gen_ele_px);
        out_tree.Branch("p4_gen_ele_py", &p4_gen_ele_py);
        out_tree.Branch("p4_gen_ele_pz", &p4_gen_ele_pz);
        out_tree.Branch("p4_gen_ele_E", &p4_gen_ele_E);
        out_tree.Branch("p4_gen_ele_motherpid", &p4_gen_ele_motherpid);
        out_tree.Branch("p4_gen_prot_px", &p4_gen_prot_px);
        out_tree.Branch("p4_gen_prot_py", &p4_gen_prot_py);
        out_tree.Branch("p4_gen_prot_pz", &p4_gen_prot_pz);
        out_tree.Branch("p4_gen_prot_E", &p4_gen_prot_E);
        out_tree.Branch("p4_gen_prot_motherpid", &p4_gen_prot_motherpid);
        out_tree.Branch("p4_gen_neutr_px", &p4_gen_neutr_px);
        out_tree.Branch("p4_gen_neutr_py", &p4_gen_neutr_py);
        out_tree.Branch("p4_gen_neutr_pz", &p4_gen_neutr_pz);
        out_tree.Branch("p4_gen_neutr_E", &p4_gen_neutr_E);
        out_tree.Branch("p4_gen_neutr_motherpid", &p4_gen_neutr_motherpid);
        out_tree.Branch("p4_gen_pip_px", &p4_gen_pip_px);
        out_tree.Branch("p4_gen_pip_py", &p4_gen_pip_py);
        out_tree.Branch("p4_gen_pip_pz", &p4_gen_pip_pz);
        out_tree.Branch("p4_gen_pip_E", &p4_gen_pip_E);
        out_tree.Branch("p4_gen_pip_motherpid", &p4_gen_pip_motherpid);
        out_tree.Branch("p4_gen_pim_px", &p4_gen_pim_px);
        out_tree.Branch("p4_gen_pim_py", &p4_gen_pim_py);
        out_tree.Branch("p4_gen_pim_pz", &p4_gen_pim_pz);
        out_tree.Branch("p4_gen_pim_E", &p4_gen_pim_E);
        out_tree.Branch("p4_gen_pim_motherpid", &p4_gen_pim_motherpid);
        out_tree.Branch("p4_gen_Kp_px", &p4_gen_Kp_px);
        out_tree.Branch("p4_gen_Kp_py", &p4_gen_Kp_py);
        out_tree.Branch("p4_gen_Kp_pz", &p4_gen_Kp_pz);
        out_tree.Branch("p4_gen_Kp_E", &p4_gen_Kp_E);
        out_tree.Branch("p4_gen_Kp_motherpid", &p4_gen_Kp_motherpid);
        out_tree.Branch("p4_gen_Km_px", &p4_gen_Km_px);
        out_tree.Branch("p4_gen_Km_py", &p4_gen_Km_py);
        out_tree.Branch("p4_gen_Km_pz", &p4_gen_Km_pz);
        out_tree.Branch("p4_gen_Km_E", &p4_gen_Km_E);
        out_tree.Branch("p4_gen_Km_motherpid", &p4_gen_Km_motherpid);
        out_tree.Branch("p4_gen_phot_px", &p4_gen_phot_px);
        out_tree.Branch("p4_gen_phot_py", &p4_gen_phot_py);
        out_tree.Branch("p4_gen_phot_pz", &p4_gen_phot_pz);
        out_tree.Branch("p4_gen_phot_E", &p4_gen_phot_E);
        out_tree.Branch("p4_gen_phot_motherpid", &p4_gen_phot_motherpid); 
        out_tree.Branch("p4_gen_pi0_px", &p4_gen_pi0_px);
        out_tree.Branch("p4_gen_pi0_py", &p4_gen_pi0_py);
        out_tree.Branch("p4_gen_pi0_pz", &p4_gen_pi0_pz);
        out_tree.Branch("p4_gen_pi0_E", &p4_gen_pi0_E);
        out_tree.Branch("p4_gen_pi0_motherpid", &p4_gen_pi0_motherpid); 
    }

    /// ///////////////////////////////////////////////////////////////
    ///  reset cut statistics
    /// ///////////////////////////////////////////////////////////////

    neg_part_count = 0;
    pos_part_count = 0;
    neut_part_count = 0;

    FD_eid_default_PID_pass = 0;
    FD_eid_charge_pass = 0;
    FD_eid_EC_outer_vs_EC_inner_pass = 0;
    FD_eid_EC_sampling_fraction_pass = 0;
    FD_eid_EC_hit_position_fiducial_pass = 0;
    FD_eid_DC_hit_position_region1_fiducial_pass = 0;
    FD_eid_DC_hit_position_region2_fiducial_pass = 0;
    FD_eid_DC_hit_position_region3_fiducial_pass = 0;
    FD_eid_DC_z_vertex_pass = 0;
    FD_eid_all_pass = 0;

    FD_protid_default_PID_pass = 0;
    FD_protid_charge_pass = 0;
    FD_protid_DC_hit_position_region1_fiducial_pass = 0;
    FD_protid_DC_hit_position_region2_fiducial_pass = 0;
    FD_protid_DC_hit_position_region3_fiducial_pass = 0;
    FD_protid_beta_pass = 0;
    FD_protid_delta_beta_pass = 0;
    FD_protid_tofmass_pass = 0;
    FD_protid_maximum_probability_pass = 0;
    FD_protid_delta_vz_pass = 0;
    FD_protid_all_pass = 0;
    FD_neutrid_default_PID_pass = 0;
    FD_neutrid_charge_pass = 0;
    FD_neutrid_beta_pass = 0;
    FD_neutrid_delta_beta_pass = 0;
    FD_neutrid_tofmass_pass = 0;
    FD_neutrid_delta_vz_pass = 0;
    FD_neutrid_all_pass = 0;
    FD_pipid_default_PID_pass = 0;
    FD_pipid_charge_pass = 0;
    FD_pipid_DC_hit_position_region1_fiducial_pass = 0;
    FD_pipid_DC_hit_position_region2_fiducial_pass = 0;
    FD_pipid_DC_hit_position_region3_fiducial_pass = 0;
    FD_pipid_beta_pass = 0;
    FD_pipid_delta_beta_pass = 0;
    FD_pipid_tofmass_pass = 0;
    FD_pipid_maximum_probability_pass = 0;
    FD_pipid_delta_vz_pass = 0;
    FD_pipid_all_pass = 0;
    FD_pimid_default_PID_pass = 0;
    FD_pimid_charge_pass = 0;
    FD_pimid_DC_hit_position_region1_fiducial_pass = 0;
    FD_pimid_DC_hit_position_region2_fiducial_pass = 0;
    FD_pimid_DC_hit_position_region3_fiducial_pass = 0;
    FD_pimid_beta_pass = 0;
    FD_pimid_delta_beta_pass = 0;
    FD_pimid_tofmass_pass = 0;
    FD_pimid_maximum_probability_pass = 0;
    FD_pimid_delta_vz_pass = 0;
    FD_pimid_all_pass = 0;
    FD_Kpid_default_PID_pass = 0;
    FD_Kpid_charge_pass = 0;
    FD_Kpid_DC_hit_position_region1_fiducial_pass = 0;
    FD_Kpid_DC_hit_position_region2_fiducial_pass = 0;
    FD_Kpid_DC_hit_position_region3_fiducial_pass = 0;
    FD_Kpid_beta_pass = 0;
    FD_Kpid_delta_beta_pass = 0;
    FD_Kpid_tofmass_pass = 0;
    FD_Kpid_maximum_probability_pass = 0;
    FD_Kpid_delta_vz_pass = 0;
    FD_Kpid_all_pass = 0;
    FD_Kmid_default_PID_pass = 0;
    FD_Kmid_charge_pass = 0;
    FD_Kmid_DC_hit_position_region1_fiducial_pass = 0;
    FD_Kmid_DC_hit_position_region2_fiducial_pass = 0;
    FD_Kmid_DC_hit_position_region3_fiducial_pass = 0;
    FD_Kmid_beta_pass = 0;
    FD_Kmid_delta_beta_pass = 0;
    FD_Kmid_tofmass_pass = 0;
    FD_Kmid_maximum_probability_pass = 0;
    FD_Kmid_delta_vz_pass = 0;
    FD_Kmid_all_pass = 0;

    CD_protid_default_PID_pass = 0;
    CD_protid_charge_pass = 0;
    CD_protid_beta_pass = 0;
    CD_protid_maximum_probability_pass = 0;
    CD_protid_delta_vz_pass = 0;
    CD_protid_all_pass = 0;
    CD_neutrid_default_PID_pass = 0;
    CD_neutrid_charge_pass = 0;
    CD_neutrid_beta_pass = 0;
    CD_neutrid_delta_vz_pass = 0;
    CD_neutrid_all_pass = 0;
    CD_pipid_default_PID_pass = 0;
    CD_pipid_charge_pass = 0;
    CD_pipid_beta_pass = 0;
    CD_pipid_maximum_probability_pass = 0;
    CD_pipid_delta_vz_pass = 0;
    CD_pipid_all_pass = 0;
    CD_pimid_default_PID_pass = 0;
    CD_pimid_charge_pass = 0;
    CD_pimid_beta_pass = 0;
    CD_pimid_maximum_probability_pass = 0;
    CD_pimid_delta_vz_pass = 0;
    CD_pimid_all_pass = 0;
    CD_Kpid_default_PID_pass = 0;
    CD_Kpid_charge_pass = 0;
    CD_Kpid_beta_pass = 0;
    CD_Kpid_maximum_probability_pass = 0;
    CD_Kpid_delta_vz_pass = 0;
    CD_Kpid_all_pass = 0;
    CD_Kmid_default_PID_pass = 0;
    CD_Kmid_charge_pass = 0;
    CD_Kmid_beta_pass = 0;
    CD_Kmid_maximum_probability_pass = 0;
    CD_Kmid_delta_vz_pass = 0;
    CD_Kmid_all_pass = 0;

    FD_photid_default_PID_pass = 0;
    FD_photid_charge_pass = 0;
    FD_photid_beta_pass = 0;
    FD_photid_EC_sampling_fraction_pass = 0;
    FD_photid_EC_hit_position_fiducial_pass = 0;
    FD_photid_all_pass = 0;

    FT_eid_charge_pass = 0;
    FT_eid_PID_pass = 0;
    FT_eid_FTCAL_fiducial_pass = 0;
    FT_eid_FTTRK_fiducial_pass = 0;
    FT_eid_FTHODO_fiducial_pass = 0;
    FT_eid_energy_vs_radius_pass = 0;
    FT_eid_all_pass = 0;

    FT_photid_charge_pass = 0;
    FT_photid_PID_pass = 0;
    FT_photid_FTCAL_fiducial_pass = 0;
    FT_photid_beta_pass = 0;
    FT_photid_all_pass = 0;

    /// ///////////////////////////////////////////////////////////////
    ///  create histograms
    /// ///////////////////////////////////////////////////////////////

    out->mkdir("event_information");
    out->cd("event_information");

    TH1F *hist_helicity;
    TH1F *hist_beam_charge;
    TH1F *hist_event_starttime;
    TH1F *hist_RF_time;
    TH1F *hist_status;

    hist_helicity = new TH1F("hist_helicity", "helicity", 5, -2.5, 2.5);
    hist_helicity->GetXaxis()->SetTitle("helicity");
    hist_helicity->GetYaxis()->SetTitle("counts");
    hist_beam_charge = new TH1F("hist_beam_charge", "beam charge", 10000, 0, 1000000);
    hist_beam_charge->GetXaxis()->SetTitle("beam charge");
    hist_beam_charge->GetYaxis()->SetTitle("counts");
    hist_event_starttime = new TH1F("hist_event_starttime", "event starttime", 600, 0, 300);
    hist_event_starttime->GetXaxis()->SetTitle("event starttime /ns");
    hist_event_starttime->GetYaxis()->SetTitle("counts");
    hist_RF_time = new TH1F("hist_RF_time", "RF time", 800, -50, 150);
    hist_RF_time->GetXaxis()->SetTitle("RF time /ns");
    hist_RF_time->GetYaxis()->SetTitle("counts");
    hist_status = new TH1F("hist_status", "status", 1800, -3000, 6000);
    hist_status->GetXaxis()->SetTitle("status");
    hist_status->GetYaxis()->SetTitle("counts");

    out->mkdir("FD_PID_electron_HTCC_plots");
    out->cd("FD_PID_electron_HTCC_plots");

    for (Int_t i = 0; i < 11; i++)
    {
        create_hist_HTCC_theta_vs_phi(i);
        create_hist_HTCC_nphe(i);
        create_hist_HTCC_nphe_vs_sampling_fraction(i);
    }

    out->mkdir("FD_PID_electron_EC_plots");
    out->cd("FD_PID_electron_EC_plots");

    for (Int_t i = 0; i < 11; i++)
    {
        create_hist_EC_PCAL_vs_EC_ECAL(i);
        create_hist_EC_outer_vs_EC_inner(i);

        create_hist_EC_total_sampling_fraction_sec1(i);
        create_hist_EC_total_sampling_fraction_sec2(i);
        create_hist_EC_total_sampling_fraction_sec3(i);
        create_hist_EC_total_sampling_fraction_sec4(i);
        create_hist_EC_total_sampling_fraction_sec5(i);
        create_hist_EC_total_sampling_fraction_sec6(i);

        create_hist_EC_PCAL_sampling_fraction_sec1(i);
        create_hist_EC_PCAL_sampling_fraction_sec2(i);
        create_hist_EC_PCAL_sampling_fraction_sec3(i);
        create_hist_EC_PCAL_sampling_fraction_sec4(i);
        create_hist_EC_PCAL_sampling_fraction_sec5(i);
        create_hist_EC_PCAL_sampling_fraction_sec6(i);

        create_hist_EC_ECAL_sampling_fraction_sec1(i);
        create_hist_EC_ECAL_sampling_fraction_sec2(i);
        create_hist_EC_ECAL_sampling_fraction_sec3(i);
        create_hist_EC_ECAL_sampling_fraction_sec4(i);
        create_hist_EC_ECAL_sampling_fraction_sec5(i);
        create_hist_EC_ECAL_sampling_fraction_sec6(i);

        create_hist_EC_PCAL_hit_position(i);
        create_hist_EC_inner_hit_position(i);
        create_hist_EC_outer_hit_position(i);
    }


    out->mkdir("FD_PID_electron_sampfrac_2D_plots");
    out->cd("FD_PID_electron_sampfrac_2D_plots");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6;

    hist_PCAL_vs_ECin_sampling_fraction_sec1 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1", "PCAL_vs_ECin_sampling_fraction_sec1", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2", "PCAL_vs_ECin_sampling_fraction_sec2", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3", "PCAL_vs_ECin_sampling_fraction_sec3", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4", "PCAL_vs_ECin_sampling_fraction_sec4", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5", "PCAL_vs_ECin_sampling_fraction_sec5", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6", "PCAL_vs_ECin_sampling_fraction_sec6", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6->GetYaxis()->SetTitle("SF_PCAL");
    
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec1;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec2;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec3;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec4;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec5;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_cut_sec6;

    hist_PCAL_vs_ECin_sampling_fraction_cut_sec1 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec1", "PCAL_vs_ECin_sampling_fraction_cut_sec1", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec1->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec1->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec2 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec2", "PCAL_vs_ECin_sampling_fraction_cut_sec2", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec2->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec2->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec3 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec3", "PCAL_vs_ECin_sampling_fraction_cut_sec3", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec3->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec3->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec4 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec4", "PCAL_vs_ECin_sampling_fraction_cut_sec4", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec4->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec4->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec5 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec5", "PCAL_vs_ECin_sampling_fraction_cut_sec5", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec5->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec5->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec6 = new TH2F("PCAL_vs_ECin_sampling_fraction_cut_sec6", "PCAL_vs_ECin_sampling_fraction_cut_sec6", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec6->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_cut_sec6->GetYaxis()->SetTitle("SF_PCAL");
    
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec1;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec2;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec3;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec4;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec5;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_pion_sec6;

    hist_PCAL_vs_ECin_sampling_fraction_pion_sec1 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec1", "PCAL_vs_ECin_sampling_fraction_pion_sec1", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec1->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec1->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec2 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec2", "PCAL_vs_ECin_sampling_fraction_pion_sec2", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec2->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec2->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec3 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec3", "PCAL_vs_ECin_sampling_fraction_pion_sec3", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec3->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec3->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec4 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec4", "PCAL_vs_ECin_sampling_fraction_pion_sec4", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec4->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec4->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec5 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec5", "PCAL_vs_ECin_sampling_fraction_pion_sec5", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec5->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec5->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec6 = new TH2F("PCAL_vs_ECin_sampling_fraction_pion_sec6", "PCAL_vs_ECin_sampling_fraction_pion_sec6", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec6->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_pion_sec6->GetYaxis()->SetTitle("SF_PCAL");
    
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec1_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec1_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p23", "PCAL_vs_ECin_sampling_fraction_sec1_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p34", "PCAL_vs_ECin_sampling_fraction_sec1_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p45", "PCAL_vs_ECin_sampling_fraction_sec1_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p56", "PCAL_vs_ECin_sampling_fraction_sec1_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p67", "PCAL_vs_ECin_sampling_fraction_sec1_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p78", "PCAL_vs_ECin_sampling_fraction_sec1_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p89", "PCAL_vs_ECin_sampling_fraction_sec1_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec1_p9", "PCAL_vs_ECin_sampling_fraction_sec1_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec1_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec2_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec2_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p23", "PCAL_vs_ECin_sampling_fraction_sec2_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p34", "PCAL_vs_ECin_sampling_fraction_sec2_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p45", "PCAL_vs_ECin_sampling_fraction_sec2_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p56", "PCAL_vs_ECin_sampling_fraction_sec2_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p67", "PCAL_vs_ECin_sampling_fraction_sec2_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p78", "PCAL_vs_ECin_sampling_fraction_sec2_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p89", "PCAL_vs_ECin_sampling_fraction_sec2_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec2_p9", "PCAL_vs_ECin_sampling_fraction_sec2_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec2_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec3_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec3_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p23", "PCAL_vs_ECin_sampling_fraction_sec3_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p34", "PCAL_vs_ECin_sampling_fraction_sec3_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p45", "PCAL_vs_ECin_sampling_fraction_sec3_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p56", "PCAL_vs_ECin_sampling_fraction_sec3_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p67", "PCAL_vs_ECin_sampling_fraction_sec3_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p78", "PCAL_vs_ECin_sampling_fraction_sec3_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p89", "PCAL_vs_ECin_sampling_fraction_sec3_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec3_p9", "PCAL_vs_ECin_sampling_fraction_sec3_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec3_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec4_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec4_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p23", "PCAL_vs_ECin_sampling_fraction_sec4_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p34", "PCAL_vs_ECin_sampling_fraction_sec4_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p45", "PCAL_vs_ECin_sampling_fraction_sec4_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p56", "PCAL_vs_ECin_sampling_fraction_sec4_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p67", "PCAL_vs_ECin_sampling_fraction_sec4_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p78", "PCAL_vs_ECin_sampling_fraction_sec4_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p89", "PCAL_vs_ECin_sampling_fraction_sec4_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec4_p9", "PCAL_vs_ECin_sampling_fraction_sec4_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec4_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec5_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec5_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p23", "PCAL_vs_ECin_sampling_fraction_sec5_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p34", "PCAL_vs_ECin_sampling_fraction_sec5_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p45", "PCAL_vs_ECin_sampling_fraction_sec5_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p56", "PCAL_vs_ECin_sampling_fraction_sec5_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p67", "PCAL_vs_ECin_sampling_fraction_sec5_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p78", "PCAL_vs_ECin_sampling_fraction_sec5_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p89", "PCAL_vs_ECin_sampling_fraction_sec5_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec5_p9", "PCAL_vs_ECin_sampling_fraction_sec5_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec5_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p23;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p34;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p45;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p56;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p67;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p78;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p89;
    TH2F *hist_PCAL_vs_ECin_sampling_fraction_sec6_p9;

    hist_PCAL_vs_ECin_sampling_fraction_sec6_p23 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p23", "PCAL_vs_ECin_sampling_fraction_sec6_p23", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p23->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p23->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p34 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p34", "PCAL_vs_ECin_sampling_fraction_sec6_p34", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p34->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p34->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p45 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p45", "PCAL_vs_ECin_sampling_fraction_sec6_p45", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p45->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p45->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p56 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p56", "PCAL_vs_ECin_sampling_fraction_sec6_p56", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p56->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p56->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p67 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p67", "PCAL_vs_ECin_sampling_fraction_sec6_p67", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p67->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p67->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p78 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p78", "PCAL_vs_ECin_sampling_fraction_sec6_p78", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p78->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p78->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p89 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p89", "PCAL_vs_ECin_sampling_fraction_sec6_p89", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p89->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p89->GetYaxis()->SetTitle("SF_PCAL");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p9 = new TH2F("PCAL_vs_ECin_sampling_fraction_sec6_p9", "PCAL_vs_ECin_sampling_fraction_sec6_p9", 200, 0, 0.2, 300, 0, 0.3);
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p9->GetXaxis()->SetTitle("SF_ECIN");
    hist_PCAL_vs_ECin_sampling_fraction_sec6_p9->GetYaxis()->SetTitle("SF_PCAL");
    
    

    out->mkdir("FD_DC_edge_plots");
    out->cd("FD_DC_edge_plots");

    TH2F *hist_DC_edge_chi2_region1_sec1_ele;
    hist_DC_edge_chi2_region1_sec1_ele = new TH2F("DC_edge_chi2_region1_sec1_ele", "DC_edge_chi2_region1_sec1_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_ele;
    hist_DC_edge_chi2_region1_sec2_ele = new TH2F("DC_edge_chi2_region1_sec2_ele", "DC_edge_chi2_region1_sec2_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_ele;
    hist_DC_edge_chi2_region1_sec3_ele = new TH2F("DC_edge_chi2_region1_sec3_ele", "DC_edge_chi2_region1_sec3_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_ele;
    hist_DC_edge_chi2_region1_sec4_ele = new TH2F("DC_edge_chi2_region1_sec4_ele", "DC_edge_chi2_region1_sec4_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_ele;
    hist_DC_edge_chi2_region1_sec5_ele = new TH2F("DC_edge_chi2_region1_sec5_ele", "DC_edge_chi2_region1_sec5_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_ele;
    hist_DC_edge_chi2_region1_sec6_ele = new TH2F("DC_edge_chi2_region1_sec6_ele", "DC_edge_chi2_region1_sec6_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_ele;
    hist_DC_edge_chi2_region2_sec1_ele = new TH2F("DC_edge_chi2_region2_sec1_ele", "DC_edge_chi2_region2_sec1_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_ele;
    hist_DC_edge_chi2_region2_sec2_ele = new TH2F("DC_edge_chi2_region2_sec2_ele", "DC_edge_chi2_region2_sec2_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_ele;
    hist_DC_edge_chi2_region2_sec3_ele = new TH2F("DC_edge_chi2_region2_sec3_ele", "DC_edge_chi2_region2_sec3_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_ele;
    hist_DC_edge_chi2_region2_sec4_ele = new TH2F("DC_edge_chi2_region2_sec4_ele", "DC_edge_chi2_region2_sec4_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_ele;
    hist_DC_edge_chi2_region2_sec5_ele = new TH2F("DC_edge_chi2_region2_sec5_ele", "DC_edge_chi2_region2_sec5_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_ele;
    hist_DC_edge_chi2_region2_sec6_ele = new TH2F("DC_edge_chi2_region2_sec6_ele", "DC_edge_chi2_region2_sec6_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_ele;
    hist_DC_edge_chi2_region3_sec1_ele = new TH2F("DC_edge_chi2_region3_sec1_ele", "DC_edge_chi2_region3_sec1_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_ele;
    hist_DC_edge_chi2_region3_sec2_ele = new TH2F("DC_edge_chi2_region3_sec2_ele", "DC_edge_chi2_region3_sec2_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_ele;
    hist_DC_edge_chi2_region3_sec3_ele = new TH2F("DC_edge_chi2_region3_sec3_ele", "DC_edge_chi2_region3_sec3_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_ele;
    hist_DC_edge_chi2_region3_sec4_ele = new TH2F("DC_edge_chi2_region3_sec4_ele", "DC_edge_chi2_region3_sec4_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_ele;
    hist_DC_edge_chi2_region3_sec5_ele = new TH2F("DC_edge_chi2_region3_sec5_ele", "DC_edge_chi2_region3_sec5_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_ele;
    hist_DC_edge_chi2_region3_sec6_ele = new TH2F("DC_edge_chi2_region3_sec6_ele", "DC_edge_chi2_region3_sec6_ele", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_ele->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_ele->GetYaxis()->SetTitle("<chi2/NDF>");
    
    
    TH2F *hist_DC_edge_chi2_region1_sec1_prot;
    hist_DC_edge_chi2_region1_sec1_prot = new TH2F("DC_edge_chi2_region1_sec1_prot", "DC_edge_chi2_region1_sec1_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_prot;
    hist_DC_edge_chi2_region1_sec2_prot = new TH2F("DC_edge_chi2_region1_sec2_prot", "DC_edge_chi2_region1_sec2_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_prot;
    hist_DC_edge_chi2_region1_sec3_prot = new TH2F("DC_edge_chi2_region1_sec3_prot", "DC_edge_chi2_region1_sec3_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_prot;
    hist_DC_edge_chi2_region1_sec4_prot = new TH2F("DC_edge_chi2_region1_sec4_prot", "DC_edge_chi2_region1_sec4_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_prot;
    hist_DC_edge_chi2_region1_sec5_prot = new TH2F("DC_edge_chi2_region1_sec5_prot", "DC_edge_chi2_region1_sec5_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_prot;
    hist_DC_edge_chi2_region1_sec6_prot = new TH2F("DC_edge_chi2_region1_sec6_prot", "DC_edge_chi2_region1_sec6_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_prot;
    hist_DC_edge_chi2_region2_sec1_prot = new TH2F("DC_edge_chi2_region2_sec1_prot", "DC_edge_chi2_region2_sec1_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_prot;
    hist_DC_edge_chi2_region2_sec2_prot = new TH2F("DC_edge_chi2_region2_sec2_prot", "DC_edge_chi2_region2_sec2_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_prot;
    hist_DC_edge_chi2_region2_sec3_prot = new TH2F("DC_edge_chi2_region2_sec3_prot", "DC_edge_chi2_region2_sec3_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_prot;
    hist_DC_edge_chi2_region2_sec4_prot = new TH2F("DC_edge_chi2_region2_sec4_prot", "DC_edge_chi2_region2_sec4_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_prot;
    hist_DC_edge_chi2_region2_sec5_prot = new TH2F("DC_edge_chi2_region2_sec5_prot", "DC_edge_chi2_region2_sec5_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_prot;
    hist_DC_edge_chi2_region2_sec6_prot = new TH2F("DC_edge_chi2_region2_sec6_prot", "DC_edge_chi2_region2_sec6_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_prot;
    hist_DC_edge_chi2_region3_sec1_prot = new TH2F("DC_edge_chi2_region3_sec1_prot", "DC_edge_chi2_region3_sec1_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_prot;
    hist_DC_edge_chi2_region3_sec2_prot = new TH2F("DC_edge_chi2_region3_sec2_prot", "DC_edge_chi2_region3_sec2_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_prot;
    hist_DC_edge_chi2_region3_sec3_prot = new TH2F("DC_edge_chi2_region3_sec3_prot", "DC_edge_chi2_region3_sec3_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_prot;
    hist_DC_edge_chi2_region3_sec4_prot = new TH2F("DC_edge_chi2_region3_sec4_prot", "DC_edge_chi2_region3_sec4_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_prot;
    hist_DC_edge_chi2_region3_sec5_prot = new TH2F("DC_edge_chi2_region3_sec5_prot", "DC_edge_chi2_region3_sec5_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_prot;
    hist_DC_edge_chi2_region3_sec6_prot = new TH2F("DC_edge_chi2_region3_sec6_prot", "DC_edge_chi2_region3_sec6_prot", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_prot->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_prot->GetYaxis()->SetTitle("<chi2/NDF>");
    
    TH2F *hist_DC_edge_chi2_region1_sec1_pip;
    hist_DC_edge_chi2_region1_sec1_pip = new TH2F("DC_edge_chi2_region1_sec1_pip", "DC_edge_chi2_region1_sec1_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_pip;
    hist_DC_edge_chi2_region1_sec2_pip = new TH2F("DC_edge_chi2_region1_sec2_pip", "DC_edge_chi2_region1_sec2_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_pip;
    hist_DC_edge_chi2_region1_sec3_pip = new TH2F("DC_edge_chi2_region1_sec3_pip", "DC_edge_chi2_region1_sec3_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_pip;
    hist_DC_edge_chi2_region1_sec4_pip = new TH2F("DC_edge_chi2_region1_sec4_pip", "DC_edge_chi2_region1_sec4_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_pip;
    hist_DC_edge_chi2_region1_sec5_pip = new TH2F("DC_edge_chi2_region1_sec5_pip", "DC_edge_chi2_region1_sec5_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_pip;
    hist_DC_edge_chi2_region1_sec6_pip = new TH2F("DC_edge_chi2_region1_sec6_pip", "DC_edge_chi2_region1_sec6_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_pip;
    hist_DC_edge_chi2_region2_sec1_pip = new TH2F("DC_edge_chi2_region2_sec1_pip", "DC_edge_chi2_region2_sec1_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_pip;
    hist_DC_edge_chi2_region2_sec2_pip = new TH2F("DC_edge_chi2_region2_sec2_pip", "DC_edge_chi2_region2_sec2_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_pip;
    hist_DC_edge_chi2_region2_sec3_pip = new TH2F("DC_edge_chi2_region2_sec3_pip", "DC_edge_chi2_region2_sec3_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_pip;
    hist_DC_edge_chi2_region2_sec4_pip = new TH2F("DC_edge_chi2_region2_sec4_pip", "DC_edge_chi2_region2_sec4_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_pip;
    hist_DC_edge_chi2_region2_sec5_pip = new TH2F("DC_edge_chi2_region2_sec5_pip", "DC_edge_chi2_region2_sec5_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_pip;
    hist_DC_edge_chi2_region2_sec6_pip = new TH2F("DC_edge_chi2_region2_sec6_pip", "DC_edge_chi2_region2_sec6_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_pip;
    hist_DC_edge_chi2_region3_sec1_pip = new TH2F("DC_edge_chi2_region3_sec1_pip", "DC_edge_chi2_region3_sec1_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_pip;
    hist_DC_edge_chi2_region3_sec2_pip = new TH2F("DC_edge_chi2_region3_sec2_pip", "DC_edge_chi2_region3_sec2_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_pip;
    hist_DC_edge_chi2_region3_sec3_pip = new TH2F("DC_edge_chi2_region3_sec3_pip", "DC_edge_chi2_region3_sec3_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_pip;
    hist_DC_edge_chi2_region3_sec4_pip = new TH2F("DC_edge_chi2_region3_sec4_pip", "DC_edge_chi2_region3_sec4_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_pip;
    hist_DC_edge_chi2_region3_sec5_pip = new TH2F("DC_edge_chi2_region3_sec5_pip", "DC_edge_chi2_region3_sec5_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_pip;
    hist_DC_edge_chi2_region3_sec6_pip = new TH2F("DC_edge_chi2_region3_sec6_pip", "DC_edge_chi2_region3_sec6_pip", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_pip->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_pip->GetYaxis()->SetTitle("<chi2/NDF>");
    
    TH2F *hist_DC_edge_chi2_region1_sec1_pim;
    hist_DC_edge_chi2_region1_sec1_pim = new TH2F("DC_edge_chi2_region1_sec1_pim", "DC_edge_chi2_region1_sec1_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_pim;
    hist_DC_edge_chi2_region1_sec2_pim = new TH2F("DC_edge_chi2_region1_sec2_pim", "DC_edge_chi2_region1_sec2_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_pim;
    hist_DC_edge_chi2_region1_sec3_pim = new TH2F("DC_edge_chi2_region1_sec3_pim", "DC_edge_chi2_region1_sec3_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_pim;
    hist_DC_edge_chi2_region1_sec4_pim = new TH2F("DC_edge_chi2_region1_sec4_pim", "DC_edge_chi2_region1_sec4_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_pim;
    hist_DC_edge_chi2_region1_sec5_pim = new TH2F("DC_edge_chi2_region1_sec5_pim", "DC_edge_chi2_region1_sec5_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_pim;
    hist_DC_edge_chi2_region1_sec6_pim = new TH2F("DC_edge_chi2_region1_sec6_pim", "DC_edge_chi2_region1_sec6_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_pim;
    hist_DC_edge_chi2_region2_sec1_pim = new TH2F("DC_edge_chi2_region2_sec1_pim", "DC_edge_chi2_region2_sec1_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_pim;
    hist_DC_edge_chi2_region2_sec2_pim = new TH2F("DC_edge_chi2_region2_sec2_pim", "DC_edge_chi2_region2_sec2_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_pim;
    hist_DC_edge_chi2_region2_sec3_pim = new TH2F("DC_edge_chi2_region2_sec3_pim", "DC_edge_chi2_region2_sec3_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_pim;
    hist_DC_edge_chi2_region2_sec4_pim = new TH2F("DC_edge_chi2_region2_sec4_pim", "DC_edge_chi2_region2_sec4_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_pim;
    hist_DC_edge_chi2_region2_sec5_pim = new TH2F("DC_edge_chi2_region2_sec5_pim", "DC_edge_chi2_region2_sec5_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_pim;
    hist_DC_edge_chi2_region2_sec6_pim = new TH2F("DC_edge_chi2_region2_sec6_pim", "DC_edge_chi2_region2_sec6_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_pim;
    hist_DC_edge_chi2_region3_sec1_pim = new TH2F("DC_edge_chi2_region3_sec1_pim", "DC_edge_chi2_region3_sec1_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_pim;
    hist_DC_edge_chi2_region3_sec2_pim = new TH2F("DC_edge_chi2_region3_sec2_pim", "DC_edge_chi2_region3_sec2_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_pim;
    hist_DC_edge_chi2_region3_sec3_pim = new TH2F("DC_edge_chi2_region3_sec3_pim", "DC_edge_chi2_region3_sec3_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_pim;
    hist_DC_edge_chi2_region3_sec4_pim = new TH2F("DC_edge_chi2_region3_sec4_pim", "DC_edge_chi2_region3_sec4_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_pim;
    hist_DC_edge_chi2_region3_sec5_pim = new TH2F("DC_edge_chi2_region3_sec5_pim", "DC_edge_chi2_region3_sec5_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_pim;
    hist_DC_edge_chi2_region3_sec6_pim = new TH2F("DC_edge_chi2_region3_sec6_pim", "DC_edge_chi2_region3_sec6_pim", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_pim->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_pim->GetYaxis()->SetTitle("<chi2/NDF>");
    
    TH2F *hist_DC_edge_chi2_region1_sec1_Kp;
    hist_DC_edge_chi2_region1_sec1_Kp = new TH2F("DC_edge_chi2_region1_sec1_Kp", "DC_edge_chi2_region1_sec1_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_Kp;
    hist_DC_edge_chi2_region1_sec2_Kp = new TH2F("DC_edge_chi2_region1_sec2_Kp", "DC_edge_chi2_region1_sec2_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_Kp;
    hist_DC_edge_chi2_region1_sec3_Kp = new TH2F("DC_edge_chi2_region1_sec3_Kp", "DC_edge_chi2_region1_sec3_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_Kp;
    hist_DC_edge_chi2_region1_sec4_Kp = new TH2F("DC_edge_chi2_region1_sec4_Kp", "DC_edge_chi2_region1_sec4_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_Kp;
    hist_DC_edge_chi2_region1_sec5_Kp = new TH2F("DC_edge_chi2_region1_sec5_Kp", "DC_edge_chi2_region1_sec5_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_Kp;
    hist_DC_edge_chi2_region1_sec6_Kp = new TH2F("DC_edge_chi2_region1_sec6_Kp", "DC_edge_chi2_region1_sec6_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_Kp;
    hist_DC_edge_chi2_region2_sec1_Kp = new TH2F("DC_edge_chi2_region2_sec1_Kp", "DC_edge_chi2_region2_sec1_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_Kp;
    hist_DC_edge_chi2_region2_sec2_Kp = new TH2F("DC_edge_chi2_region2_sec2_Kp", "DC_edge_chi2_region2_sec2_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_Kp;
    hist_DC_edge_chi2_region2_sec3_Kp = new TH2F("DC_edge_chi2_region2_sec3_Kp", "DC_edge_chi2_region2_sec3_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_Kp;
    hist_DC_edge_chi2_region2_sec4_Kp = new TH2F("DC_edge_chi2_region2_sec4_Kp", "DC_edge_chi2_region2_sec4_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_Kp;
    hist_DC_edge_chi2_region2_sec5_Kp = new TH2F("DC_edge_chi2_region2_sec5_Kp", "DC_edge_chi2_region2_sec5_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_Kp;
    hist_DC_edge_chi2_region2_sec6_Kp = new TH2F("DC_edge_chi2_region2_sec6_Kp", "DC_edge_chi2_region2_sec6_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_Kp;
    hist_DC_edge_chi2_region3_sec1_Kp = new TH2F("DC_edge_chi2_region3_sec1_Kp", "DC_edge_chi2_region3_sec1_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_Kp;
    hist_DC_edge_chi2_region3_sec2_Kp = new TH2F("DC_edge_chi2_region3_sec2_Kp", "DC_edge_chi2_region3_sec2_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_Kp;
    hist_DC_edge_chi2_region3_sec3_Kp = new TH2F("DC_edge_chi2_region3_sec3_Kp", "DC_edge_chi2_region3_sec3_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_Kp;
    hist_DC_edge_chi2_region3_sec4_Kp = new TH2F("DC_edge_chi2_region3_sec4_Kp", "DC_edge_chi2_region3_sec4_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_Kp;
    hist_DC_edge_chi2_region3_sec5_Kp = new TH2F("DC_edge_chi2_region3_sec5_Kp", "DC_edge_chi2_region3_sec5_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_Kp;
    hist_DC_edge_chi2_region3_sec6_Kp = new TH2F("DC_edge_chi2_region3_sec6_Kp", "DC_edge_chi2_region3_sec6_Kp", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_Kp->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_Kp->GetYaxis()->SetTitle("<chi2/NDF>");
    
    TH2F *hist_DC_edge_chi2_region1_sec1_Km;
    hist_DC_edge_chi2_region1_sec1_Km = new TH2F("DC_edge_chi2_region1_sec1_Km", "DC_edge_chi2_region1_sec1_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec1_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec1_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec2_Km;
    hist_DC_edge_chi2_region1_sec2_Km = new TH2F("DC_edge_chi2_region1_sec2_Km", "DC_edge_chi2_region1_sec2_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec2_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec2_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec3_Km;
    hist_DC_edge_chi2_region1_sec3_Km = new TH2F("DC_edge_chi2_region1_sec3_Km", "DC_edge_chi2_region1_sec3_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec3_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec3_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec4_Km;
    hist_DC_edge_chi2_region1_sec4_Km = new TH2F("DC_edge_chi2_region1_sec4_Km", "DC_edge_chi2_region1_sec4_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec4_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec4_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec5_Km;
    hist_DC_edge_chi2_region1_sec5_Km = new TH2F("DC_edge_chi2_region1_sec5_Km", "DC_edge_chi2_region1_sec5_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec5_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec5_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region1_sec6_Km;
    hist_DC_edge_chi2_region1_sec6_Km = new TH2F("DC_edge_chi2_region1_sec6_Km", "DC_edge_chi2_region1_sec6_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region1_sec6_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region1_sec6_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec1_Km;
    hist_DC_edge_chi2_region2_sec1_Km = new TH2F("DC_edge_chi2_region2_sec1_Km", "DC_edge_chi2_region2_sec1_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec1_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec1_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec2_Km;
    hist_DC_edge_chi2_region2_sec2_Km = new TH2F("DC_edge_chi2_region2_sec2_Km", "DC_edge_chi2_region2_sec2_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec2_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec2_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec3_Km;
    hist_DC_edge_chi2_region2_sec3_Km = new TH2F("DC_edge_chi2_region2_sec3_Km", "DC_edge_chi2_region2_sec3_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec3_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec3_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec4_Km;
    hist_DC_edge_chi2_region2_sec4_Km = new TH2F("DC_edge_chi2_region2_sec4_Km", "DC_edge_chi2_region2_sec4_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec4_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec4_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec5_Km;
    hist_DC_edge_chi2_region2_sec5_Km = new TH2F("DC_edge_chi2_region2_sec5_Km", "DC_edge_chi2_region2_sec5_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec5_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec5_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region2_sec6_Km;
    hist_DC_edge_chi2_region2_sec6_Km = new TH2F("DC_edge_chi2_region2_sec6_Km", "DC_edge_chi2_region2_sec6_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region2_sec6_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region2_sec6_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec1_Km;
    hist_DC_edge_chi2_region3_sec1_Km = new TH2F("DC_edge_chi2_region3_sec1_Km", "DC_edge_chi2_region3_sec1_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec1_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec1_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec2_Km;
    hist_DC_edge_chi2_region3_sec2_Km = new TH2F("DC_edge_chi2_region3_sec2_Km", "DC_edge_chi2_region3_sec2_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec2_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec2_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec3_Km;
    hist_DC_edge_chi2_region3_sec3_Km = new TH2F("DC_edge_chi2_region3_sec3_Km", "DC_edge_chi2_region3_sec3_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec3_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec3_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec4_Km;
    hist_DC_edge_chi2_region3_sec4_Km = new TH2F("DC_edge_chi2_region3_sec4_Km", "DC_edge_chi2_region3_sec4_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec4_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec4_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec5_Km;
    hist_DC_edge_chi2_region3_sec5_Km = new TH2F("DC_edge_chi2_region3_sec5_Km", "DC_edge_chi2_region3_sec5_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec5_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec5_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    TH2F *hist_DC_edge_chi2_region3_sec6_Km;
    hist_DC_edge_chi2_region3_sec6_Km = new TH2F("DC_edge_chi2_region3_sec6_Km", "DC_edge_chi2_region3_sec6_Km", 300, 0, 30, 5000, 0, 500);
    hist_DC_edge_chi2_region3_sec6_Km->GetXaxis()->SetTitle("<edge>");
    hist_DC_edge_chi2_region3_sec6_Km->GetYaxis()->SetTitle("<chi2/NDF>");
    
    


    out->mkdir("FD_PID_electron_DC_plots");
    out->cd("FD_PID_electron_DC_plots");

    for (Int_t i = 0; i < 11; i++)
    {

        create_hist_DC_hit_position_region1(i);
        create_hist_DC_hit_position_region2(i);
        create_hist_DC_hit_position_region3(i);

        create_hist_DC_z_vertex_sec1(i);
        create_hist_DC_z_vertex_sec2(i);
        create_hist_DC_z_vertex_sec3(i);
        create_hist_DC_z_vertex_sec4(i);
        create_hist_DC_z_vertex_sec5(i);
        create_hist_DC_z_vertex_sec6(i);
    }

    TH2F *hist_DC_hit_position_region2_cut5a;
    hist_DC_hit_position_region2_cut5a = new TH2F("DC_hit_position_region2_cut_05a", "DC_hit_position_region2_cut_05a", 1000, -500, 500, 1000, -500, 500);
    hist_DC_hit_position_region2_cut5a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_cut5a->GetYaxis()->SetTitle("y /cm");

    out->mkdir("FD_PID_hadron_DC_fiducial_plot");
    out->cd("FD_PID_hadron_DC_fiducial_plot");

    for (Int_t i = 0; i < 60; i++)
    {
        create_DC_hit_position_region1_hadron(i);
    }
    for (Int_t i = 0; i < 60; i++)
    {
        create_DC_hit_position_region2_hadron(i);
    }
    for (Int_t i = 0; i < 60; i++)
    {
        create_DC_hit_position_region3_hadron(i);
    }

    TH2F *hist_DC_hit_position_region2_hadron_cut_02a;
    TH2F *hist_DC_hit_position_region2_hadron_cut_12a;
    TH2F *hist_DC_hit_position_region2_hadron_cut_22a;
    TH2F *hist_DC_hit_position_region2_hadron_cut_32a;
    TH2F *hist_DC_hit_position_region2_hadron_cut_42a;
    TH2F *hist_DC_hit_position_region2_hadron_cut_52a;

    hist_DC_hit_position_region2_hadron_cut_02a = new TH2F("DC_hit_position_region2_hadron_cut_02a", "DC_hit_position_region2_hadron_cut_02a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_02a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_02a->GetYaxis()->SetTitle("y /cm");
    hist_DC_hit_position_region2_hadron_cut_12a = new TH2F("DC_hit_position_region2_hadron_cut_12a", "DC_hit_position_region2_hadron_cut_12a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_12a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_12a->GetYaxis()->SetTitle("y /cm");
    hist_DC_hit_position_region2_hadron_cut_22a = new TH2F("DC_hit_position_region2_hadron_cut_22a", "DC_hit_position_region2_hadron_cut_22a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_22a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_22a->GetYaxis()->SetTitle("y /cm");
    hist_DC_hit_position_region2_hadron_cut_32a = new TH2F("DC_hit_position_region2_hadron_cut_32a", "DC_hit_position_region2_hadron_cut_32a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_32a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_32a->GetYaxis()->SetTitle("y /cm");
    hist_DC_hit_position_region2_hadron_cut_42a = new TH2F("DC_hit_position_region2_hadron_cut_42a", "DC_hit_position_region2_hadron_cut_42a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_42a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_42a->GetYaxis()->SetTitle("y /cm");
    hist_DC_hit_position_region2_hadron_cut_52a = new TH2F("DC_hit_position_region2_hadron_cut_52a", "DC_hit_position_region2_hadron_cut_52a", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region2_hadron_cut_52a->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron_cut_52a->GetYaxis()->SetTitle("y /cm");

    out->mkdir("FD_PID_hadron_EC_plots");
    out->cd("FD_PID_hadron_EC_plots");

    for (Int_t i = 0; i < 60; i++)
    {
        create_hist_EC_outer_vs_EC_inner_hadron(i);
    }

    out->mkdir("FD_PID_hadron_beta_plots");
    out->cd("FD_PID_hadron_beta_plots");

    for (Int_t i = 0; i < 60; i++)
    {
        create_hist_beta_vs_p(i);
        create_hist_beta_vs_p_sec1(i);
        create_hist_beta_vs_p_sec2(i);
        create_hist_beta_vs_p_sec3(i);
        create_hist_beta_vs_p_sec4(i);
        create_hist_beta_vs_p_sec5(i);
        create_hist_beta_vs_p_sec6(i);
    }

    for (Int_t i = 0; i < 60; i++)
    {
        create_hist_delta_beta_vs_p(i);
    }

    for (Int_t i = 0; i < 60; i++)
    {
        create_hist_delta_beta(i);
    }


    out->mkdir("FD_PID_hadron_vertex_plots");
    out->cd("FD_PID_hadron_vertex_plots");

    for (Int_t i = 0; i < 60; i++)
    {
        create_hist_delta_vz(i);
    }

    // CD hadrons

    out->mkdir("CD_PID_hadron_beta_plots");
    out->cd("CD_PID_hadron_beta_plots");

    for (Int_t i = 0; i < 6; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }
    for (Int_t i = 10; i < 16; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }
    for (Int_t i = 20; i < 26; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }
    for (Int_t i = 30; i < 36; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }
    for (Int_t i = 40; i < 46; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }
    for (Int_t i = 50; i < 56; i++)
    {
        create_hist_CD_beta_vs_p(i);
    }

    out->mkdir("CD_PID_hadron_vertex_plots");
    out->cd("CD_PID_hadron_vertex_plots");

    for (Int_t i = 0; i < 6; i++)
    {
        create_hist_CD_delta_vz(i);
    }
    for (Int_t i = 10; i < 16; i++)
    {
        create_hist_CD_delta_vz(i);
    }
    for (Int_t i = 20; i < 26; i++)
    {
        create_hist_CD_delta_vz(i);
    }
    for (Int_t i = 30; i < 36; i++)
    {
        create_hist_CD_delta_vz(i);
    }
    for (Int_t i = 40; i < 46; i++)
    {
        create_hist_CD_delta_vz(i);
    }
    for (Int_t i = 50; i < 56; i++)
    {
        create_hist_CD_delta_vz(i);
    }

    out->mkdir("FD_PID_neutrals_plots");
    out->cd("FD_PID_neutrals_plots");

    for (Int_t i = 0; i <= 4; i++)
    {
        create_hist_beta_vs_p_phot(i);
        create_hist_beta_phot(i);
        create_hist_EC_sampling_fraction_phot(i);
        create_hist_EC_PCAL_vs_EC_ECAL_phot(i);
        create_hist_EC_PCAL_hit_position_phot(i);
    }

    out->mkdir("FT_PID_plots");
    out->cd("FT_PID_plots");

    // elctrons
    for (Int_t i = 0; i < 8; i++)
    {
        create_hist_FT_FTCAL_energy_vs_radius(i);
        create_hist_FT_FTCAL_hit_position(i);
        create_hist_FT_FTTRK_hit_position(i);
        create_hist_FT_FTHODO_hit_position(i);
        create_hist_FT_beta(i);
    }

    // photons
    for (Int_t i = 10; i < 16; i++)
    {
        create_hist_FT_FTCAL_energy_vs_radius(i);
        create_hist_FT_FTCAL_hit_position(i);
        create_hist_FT_FTTRK_hit_position(i);
        create_hist_FT_FTHODO_hit_position(i);
        create_hist_FT_beta(i);
    }



    /// //////////////////////////////////////////////////////////////////////////////
    /// Fiducial cut plots

    out->mkdir("EC_fiducial_cut_plots");
    out->cd("EC_fiducial_cut_plots");

    TH2F *hist_electron_sampfrac_vs_u_coord_sec1;
    TH2F *hist_electron_sampfrac_vs_u_coord_sec2;
    TH2F *hist_electron_sampfrac_vs_u_coord_sec3;
    TH2F *hist_electron_sampfrac_vs_u_coord_sec4;
    TH2F *hist_electron_sampfrac_vs_u_coord_sec5;
    TH2F *hist_electron_sampfrac_vs_u_coord_sec6;
    TH2F *hist_electron_sampfrac_vs_u_coord;

    TH2F *hist_electron_sampfrac_vs_v_coord_sec1;
    TH2F *hist_electron_sampfrac_vs_v_coord_sec2;
    TH2F *hist_electron_sampfrac_vs_v_coord_sec3;
    TH2F *hist_electron_sampfrac_vs_v_coord_sec4;
    TH2F *hist_electron_sampfrac_vs_v_coord_sec5;
    TH2F *hist_electron_sampfrac_vs_v_coord_sec6;
    TH2F *hist_electron_sampfrac_vs_v_coord;

    TH2F *hist_electron_sampfrac_vs_w_coord_sec1;
    TH2F *hist_electron_sampfrac_vs_w_coord_sec2;
    TH2F *hist_electron_sampfrac_vs_w_coord_sec3;
    TH2F *hist_electron_sampfrac_vs_w_coord_sec4;
    TH2F *hist_electron_sampfrac_vs_w_coord_sec5;
    TH2F *hist_electron_sampfrac_vs_w_coord_sec6;
    TH2F *hist_electron_sampfrac_vs_w_coord;

    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec1;
    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec2;
    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec3;
    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec4;
    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec5;
    TH3F *hist_electron_sampfrac_vs_u_coord_vs_p_sec6;

    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec1;
    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec2;
    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec3;
    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec4;
    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec5;
    TH3F *hist_electron_sampfrac_vs_v_coord_vs_p_sec6;

    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec1;
    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec2;
    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec3;
    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec4;
    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec5;
    TH3F *hist_electron_sampfrac_vs_w_coord_vs_p_sec6;

    hist_electron_sampfrac_vs_u_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec1", "electron E/p vs u coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec2", "electron E/p vs u coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec3", "electron E/p vs u coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec4", "electron E/p vs u coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec5", "electron E/p vs u coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_u_coord_vs_p_sec6", "electron E/p vs u coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");

    hist_electron_sampfrac_vs_v_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec1", "electron E/p vs v coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec2", "electron E/p vs v coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec3", "electron E/p vs v coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec4", "electron E/p vs v coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec5", "electron E/p vs v coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_v_coord_vs_p_sec6", "electron E/p vs v coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");

    hist_electron_sampfrac_vs_w_coord_vs_p_sec1 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec1", "electron E/p vs w coordinate vs p sector 1", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec1->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec2 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec2", "electron E/p vs w coordinate vs p sector 2", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec2->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec3 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec3", "electron E/p vs w coordinate vs p sector 3", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec3->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec4 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec4", "electron E/p vs w coordinate vs p sector 4", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec4->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec5 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec5", "electron E/p vs w coordinate vs p sector 5", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec5->GetZaxis()->SetTitle("p /GeV");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec6 = new TH3F("hist_electron_sampfrac_vs_w_coord_vs_p_sec6", "electron E/p vs w coordinate vs p sector 6", 450, 0, 450, 50, 0, 0.5, 100, 0.0, Ebeam + 1);
    hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_vs_p_sec6->GetZaxis()->SetTitle("p /GeV");

    hist_electron_sampfrac_vs_u_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec1", "electron E/p vs u coordinate sector 1", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec1->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec2", "electron E/p vs u coordinate sector 2", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec2->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec3", "electron E/p vs u coordinate sector 3", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec3->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec4", "electron E/p vs u coordinate sector 4", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec4->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec5", "electron E/p vs u coordinate sector 5", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec5->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_u_coord_sec6", "electron E/p vs u coordinate sector 6", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord_sec6->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_u_coord = new TH2F("hist_electron_sampfrac_vs_u_coord", "electron E/p vs u coordinate", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_u_coord->GetXaxis()->SetTitle("u /cm");
    hist_electron_sampfrac_vs_u_coord->GetYaxis()->SetTitle("E/p");

    hist_electron_sampfrac_vs_v_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec1", "electron E/p vs v coordinate sector 1", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec1->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec2", "electron E/p vs v coordinate sector 2", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec2->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec3", "electron E/p vs v coordinate sector 3", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec3->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec4", "electron E/p vs v coordinate sector 4", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec4->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec5", "electron E/p vs v coordinate sector 5", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec5->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_v_coord_sec6", "electron E/p vs v coordinate sector 6", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord_sec6->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_v_coord = new TH2F("hist_electron_sampfrac_vs_v_coord", "electron E/p vs v coordinate", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_v_coord->GetXaxis()->SetTitle("v /cm");
    hist_electron_sampfrac_vs_v_coord->GetYaxis()->SetTitle("E/p");

    hist_electron_sampfrac_vs_w_coord_sec1 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec1", "electron E/p vs w coordinate sector 1", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec1->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec1->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_sec2 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec2", "electron E/p vs w coordinate sector 2", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec2->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec2->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_sec3 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec3", "electron E/p vs w coordinate sector 3", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec3->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec3->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_sec4 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec4", "electron E/p vs w coordinate sector 4", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec4->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec4->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_sec5 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec5", "electron E/p vs w coordinate sector 5", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec5->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec5->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord_sec6 = new TH2F("hist_electron_sampfrac_vs_w_coord_sec6", "electron E/p vs w coordinate sector 6", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord_sec6->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord_sec6->GetYaxis()->SetTitle("E/p");
    hist_electron_sampfrac_vs_w_coord = new TH2F("hist_electron_sampfrac_vs_w_coord", "electron E/p vs w coordinate", 900, 0, 450, 50, 0, 0.5);
    hist_electron_sampfrac_vs_w_coord->GetXaxis()->SetTitle("w /cm");
    hist_electron_sampfrac_vs_w_coord->GetYaxis()->SetTitle("E/p");

    out->mkdir("EC_fiducial_cut_photons");
    out->cd("EC_fiducial_cut_photons");

    TH1F *hist_photon_u_coord_sec1;
    TH1F *hist_photon_u_coord_sec2;
    TH1F *hist_photon_u_coord_sec3;
    TH1F *hist_photon_u_coord_sec4;
    TH1F *hist_photon_u_coord_sec5;
    TH1F *hist_photon_u_coord_sec6;
    TH1F *hist_photon_u_coord;

    TH1F *hist_photon_v_coord_sec1;
    TH1F *hist_photon_v_coord_sec2;
    TH1F *hist_photon_v_coord_sec3;
    TH1F *hist_photon_v_coord_sec4;
    TH1F *hist_photon_v_coord_sec5;
    TH1F *hist_photon_v_coord_sec6;
    TH1F *hist_photon_v_coord;

    TH1F *hist_photon_w_coord_sec1;
    TH1F *hist_photon_w_coord_sec2;
    TH1F *hist_photon_w_coord_sec3;
    TH1F *hist_photon_w_coord_sec4;
    TH1F *hist_photon_w_coord_sec5;
    TH1F *hist_photon_w_coord_sec6;
    TH1F *hist_photon_w_coord;

    hist_photon_u_coord_sec1 = new TH1F("hist_photon_u_coord_sec1", "photon u coordinate sector 1", 900, 0, 450);
    hist_photon_u_coord_sec2 = new TH1F("hist_photon_u_coord_sec2", "photon u coordinate sector 2", 900, 0, 450);
    hist_photon_u_coord_sec3 = new TH1F("hist_photon_u_coord_sec3", "photon u coordinate sector 3", 900, 0, 450);
    hist_photon_u_coord_sec4 = new TH1F("hist_photon_u_coord_sec4", "photon u coordinate sector 3", 900, 0, 450);
    hist_photon_u_coord_sec5 = new TH1F("hist_photon_u_coord_sec5", "photon u coordinate sector 4", 900, 0, 450);
    hist_photon_u_coord_sec6 = new TH1F("hist_photon_u_coord_sec6", "photon u coordinate sector 5", 900, 0, 450);
    hist_photon_u_coord = new TH1F("hist_photon_u_coord", "photon u coordinate", 900, 0, 450);

    hist_photon_v_coord_sec1 = new TH1F("hist_photon_v_coord_sec1", "photon v coordinate sector 1", 900, 0, 450);
    hist_photon_v_coord_sec2 = new TH1F("hist_photon_v_coord_sec2", "photon v coordinate sector 2", 900, 0, 450);
    hist_photon_v_coord_sec3 = new TH1F("hist_photon_v_coord_sec3", "photon v coordinate sector 3", 900, 0, 450);
    hist_photon_v_coord_sec4 = new TH1F("hist_photon_v_coord_sec4", "photon v coordinate sector 3", 900, 0, 450);
    hist_photon_v_coord_sec5 = new TH1F("hist_photon_v_coord_sec5", "photon v coordinate sector 4", 900, 0, 450);
    hist_photon_v_coord_sec6 = new TH1F("hist_photon_v_coord_sec6", "photon v coordinate sector 5", 900, 0, 450);
    hist_photon_v_coord = new TH1F("hist_photon_v_coord", "photon v coordinate", 900, 0, 450);

    hist_photon_w_coord_sec1 = new TH1F("hist_photon_w_coord_sec1", "photon w coordinate sector 1", 900, 0, 450);
    hist_photon_w_coord_sec2 = new TH1F("hist_photon_w_coord_sec2", "photon w coordinate sector 2", 900, 0, 450);
    hist_photon_w_coord_sec3 = new TH1F("hist_photon_w_coord_sec3", "photon w coordinate sector 3", 900, 0, 450);
    hist_photon_w_coord_sec4 = new TH1F("hist_photon_w_coord_sec4", "photon w coordinate sector 3", 900, 0, 450);
    hist_photon_w_coord_sec5 = new TH1F("hist_photon_w_coord_sec5", "photon w coordinate sector 4", 900, 0, 450);
    hist_photon_w_coord_sec6 = new TH1F("hist_photon_w_coord_sec6", "photon w coordinate sector 5", 900, 0, 450);
    hist_photon_w_coord = new TH1F("hist_photon_w_coord", "photon w coordinate", 900, 0, 450);

    /// //////////////////////////////////////////////////////////////////////////////
    /// new fiducial cut plots

    out->mkdir("fiducial_cuts_new");
    out->cd("fiducial_cuts_new");

    TH2F *hist_electron_sampfrac_vs_phi_sec1[30];
    TH2F *hist_electron_sampfrac_vs_phi_sec2[30];
    TH2F *hist_electron_sampfrac_vs_phi_sec3[30];
    TH2F *hist_electron_sampfrac_vs_phi_sec4[30];
    TH2F *hist_electron_sampfrac_vs_phi_sec5[30];
    TH2F *hist_electron_sampfrac_vs_phi_sec6[30];
    TH2F *hist_electron_sampfrac_vs_phi[30];

    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec1_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec1_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec1[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec1[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec1[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec2_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec2_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec2[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec2[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec2[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec3_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec3_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec3[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec3[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec3[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec4_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec4_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec4[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec4[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec4[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec5_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec5_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec5[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec5[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec5[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_sec6_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_sec6_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi_sec6[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi_sec6[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi_sec6[j]->GetYaxis()->SetTitle("E/p");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_electron_sampfrac_vs_phi_thetabin_%01d", j);
        sprintf(title, "hist_electron_sampfrac_vs_phi_thetabin_%01d", j);
        hist_electron_sampfrac_vs_phi[j] = new TH2F(name, title, 400, -40, 40, 50, 0, 0.5);
        hist_electron_sampfrac_vs_phi[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_electron_sampfrac_vs_phi[j]->GetYaxis()->SetTitle("E/p");
    }

    TH1F *hist_DC_region1_electron_counts_vs_phi_sec1[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi_sec2[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi_sec3[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi_sec4[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi_sec5[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi_sec6[30];
    TH1F *hist_DC_region1_electron_counts_vs_phi[30];

    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec1[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec1[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec1[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec2[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec2[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec2[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec3[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec3[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec3[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec4[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec4[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec4[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec5[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec5[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec5[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi_sec6[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi_sec6[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi_sec6[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region1_electron_counts_vs_phi_thetabin_%01d", j);
        sprintf(title, "hist_DC_region1_electron_counts_vs_phi_thetabin_%01d", j);
        hist_DC_region1_electron_counts_vs_phi[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region1_electron_counts_vs_phi[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region1_electron_counts_vs_phi[j]->GetYaxis()->SetTitle("counts");
    }

    TH1F *hist_DC_region2_electron_counts_vs_phi_sec1[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi_sec2[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi_sec3[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi_sec4[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi_sec5[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi_sec6[30];
    TH1F *hist_DC_region2_electron_counts_vs_phi[30];

    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec1[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec1[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec1[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec2[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec2[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec2[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec3[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec3[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec3[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec4[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec4[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec4[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec5[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec5[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec5[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi_sec6[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi_sec6[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi_sec6[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region2_electron_counts_vs_phi_thetabin_%01d", j);
        sprintf(title, "hist_DC_region2_electron_counts_vs_phi_thetabin_%01d", j);
        hist_DC_region2_electron_counts_vs_phi[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region2_electron_counts_vs_phi[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region2_electron_counts_vs_phi[j]->GetYaxis()->SetTitle("counts");
    }

    TH1F *hist_DC_region3_electron_counts_vs_phi_sec1[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi_sec2[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi_sec3[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi_sec4[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi_sec5[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi_sec6[30];
    TH1F *hist_DC_region3_electron_counts_vs_phi[30];

    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec1_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec1[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec1[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec1[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec2_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec2[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec2[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec2[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec3_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec3[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec3[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec3[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec4_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec4[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec4[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec4[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec5_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec5[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec5[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec5[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_sec6_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi_sec6[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi_sec6[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi_sec6[j]->GetYaxis()->SetTitle("counts");
    }
    for (Int_t j = 0; j < 30; j++)
    {
        sprintf(name, "hist_DC_region3_electron_counts_vs_phi_thetabin_%01d", j);
        sprintf(title, "hist_DC_region3_electron_counts_vs_phi_thetabin_%01d", j);
        hist_DC_region3_electron_counts_vs_phi[j] = new TH1F(name, title, 400, -40, 40);
        hist_DC_region3_electron_counts_vs_phi[j]->GetXaxis()->SetTitle("#phi_{local}");
        hist_DC_region3_electron_counts_vs_phi[j]->GetYaxis()->SetTitle("counts");
    }

    TH2F *hist_DC_hit_position_region1_new_cut;
    TH2F *hist_DC_hit_position_region2_new_cut;
    TH2F *hist_DC_hit_position_region3_new_cut;

    hist_DC_hit_position_region1_new_cut = new TH2F("DC_hit_position_region1_new_cut", "DC_hit_position_region1_new_cut", 1000, -450, 450, 1000, -450, 450);
    hist_DC_hit_position_region1_new_cut->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region1_new_cut->GetYaxis()->SetTitle("y /cm");

    hist_DC_hit_position_region2_new_cut = new TH2F("DC_hit_position_region2_new_cut", "DC_hit_position_region2_new_cut", 1000, -500, 500, 1000, -500, 500);
    hist_DC_hit_position_region2_new_cut->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_new_cut->GetYaxis()->SetTitle("y /cm");

    hist_DC_hit_position_region3_new_cut = new TH2F("DC_hit_position_region3_new_cut", "DC_hit_position_region3_new_cut", 1000, -500, 500, 1000, -500, 500);
    hist_DC_hit_position_region3_new_cut->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region3_new_cut->GetYaxis()->SetTitle("y /cm");

    /// //////////////////////////////////////////////////////////////////////////////
    /// Fiducial cut plots

    out->mkdir("HTCC_Nphe");
    out->cd("HTCC_Nphe");

    TH2F *hist_HTCC_Nphe_vs_momentum;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec1;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec2;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec3;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec4;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec5;
    TH2F *hist_HTCC_Nphe_vs_momentum_sec6;

    TH2F *hist_HTCC_Nphe_vs_momentum_prot;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec1;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec2;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec3;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec4;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec5;
    TH2F *hist_HTCC_Nphe_vs_momentum_prot_sec6;

    TH2F *hist_HTCC_Nphe_vs_momentum_pip;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec1;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec2;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec3;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec4;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec5;
    TH2F *hist_HTCC_Nphe_vs_momentum_pip_sec6;

    TH2F *hist_HTCC_Nphe_vs_beta;
    TH2F *hist_HTCC_Nphe_vs_beta_sec1;
    TH2F *hist_HTCC_Nphe_vs_beta_sec2;
    TH2F *hist_HTCC_Nphe_vs_beta_sec3;
    TH2F *hist_HTCC_Nphe_vs_beta_sec4;
    TH2F *hist_HTCC_Nphe_vs_beta_sec5;
    TH2F *hist_HTCC_Nphe_vs_beta_sec6;

    TH2F *hist_HTCC_Nphe_vs_beta_prot;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec1;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec2;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec3;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec4;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec5;
    TH2F *hist_HTCC_Nphe_vs_beta_prot_sec6;

    TH2F *hist_HTCC_Nphe_vs_beta_pip;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec1;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec2;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec3;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec4;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec5;
    TH2F *hist_HTCC_Nphe_vs_beta_pip_sec6;

    hist_HTCC_Nphe_vs_momentum = new TH2F("hist_HTCC_Nphe_vs_momentum", "HTCC Nphe vs momentum", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec1", "HTCC Nphe vs momentum sector 1", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec1->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec1->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec2", "HTCC Nphe vs momentum sector 2", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec2->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec2->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec3", "HTCC Nphe vs momentum sector 3", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec3->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec3->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec4", "HTCC Nphe vs momentum sector 4", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec4->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec4->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec5", "HTCC Nphe vs momentum sector 5", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec5->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec5->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_sec6", "HTCC Nphe vs momentum sector 6", 500, 0, Ebeam + 0.5, 200, 0, 200);
    hist_HTCC_Nphe_vs_momentum_sec6->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_sec6->GetYaxis()->SetTitle("Nphe");

    hist_HTCC_Nphe_vs_momentum_prot = new TH2F("hist_HTCC_Nphe_vs_momentum_prot", "HTCC Nphe vs momentum for protons", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec1", "HTCC Nphe vs momentum for protons sector 1", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec1->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec1->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec2", "HTCC Nphe vs momentum for protons sector 2", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec2->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec2->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec3", "HTCC Nphe vs momentum for protons sector 3", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec3->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec3->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec4", "HTCC Nphe vs momentum for protons sector 4", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec4->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec4->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec5", "HTCC Nphe vs momentum for protons sector 5", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec5->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec5->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_prot_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_prot_sec6", "HTCC Nphe vs momentum for protons sector 6", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_prot_sec6->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_prot_sec6->GetYaxis()->SetTitle("Nphe");

    hist_HTCC_Nphe_vs_momentum_pip = new TH2F("hist_HTCC_Nphe_vs_momentum_pip", "HTCC Nphe vs momentum for pip", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec1 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec1", "HTCC Nphe vs momentum for pip sector 1", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec1->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec1->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec2 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec2", "HTCC Nphe vs momentum for pip sector 2", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec2->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec2->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec3 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec3", "HTCC Nphe vs momentum for pip sector 3", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec3->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec3->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec4 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec4", "HTCC Nphe vs momentum for pip sector 4", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec4->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec4->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec5 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec5", "HTCC Nphe vs momentum for pip sector 5", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec5->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec5->GetYaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_momentum_pip_sec6 = new TH2F("hist_HTCC_Nphe_vs_momentum_pip_sec6", "HTCC Nphe vs momentum for pip sector 6", 500, 0, Ebeam - 0.5, 60, 0, 60);
    hist_HTCC_Nphe_vs_momentum_pip_sec6->GetXaxis()->SetTitle("p /GeV");
    hist_HTCC_Nphe_vs_momentum_pip_sec6->GetYaxis()->SetTitle("Nphe");

    hist_HTCC_Nphe_vs_beta = new TH2F("hist_HTCC_Nphe_vs_beta", "HTCC Nphe vs beta", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_sec1", "HTCC Nphe vs beta sector 1", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec1->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec1->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_sec2", "HTCC Nphe vs beta sector 2", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec2->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec2->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_sec3", "HTCC Nphe vs beta sector 3", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec3->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec3->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_sec4", "HTCC Nphe vs beta sector 4", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec4->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec4->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_sec5", "HTCC Nphe vs beta sector 5", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec5->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec5->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_sec6", "HTCC Nphe vs beta sector 6", 300, 0, 300, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_sec6->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_sec6->GetYaxis()->SetTitle("#beta");

    hist_HTCC_Nphe_vs_beta_prot = new TH2F("hist_HTCC_Nphe_vs_beta_prot", "HTCC Nphe vs beta for protons", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec1", "HTCC Nphe vs beta for protons sector 1", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec1->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec1->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec2", "HTCC Nphe vs beta for protons sector 2", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec2->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec2->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec3", "HTCC Nphe vs beta for protons sector 3", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec3->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec3->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec4", "HTCC Nphe vs beta for protons sector 4", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec4->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec4->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec5", "HTCC Nphe vs beta for protons sector 5", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec5->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec5->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_prot_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_prot_sec6", "HTCC Nphe vs beta for protons sector 6", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_prot_sec6->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_prot_sec6->GetYaxis()->SetTitle("#beta");

    hist_HTCC_Nphe_vs_beta_pip = new TH2F("hist_HTCC_Nphe_vs_beta_pip", "HTCC Nphe vs beta for pip", 100, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec1 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec1", "HTCC Nphe vs beta for pip sector 1", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec1->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec1->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec2 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec2", "HTCC Nphe vs beta for pip sector 2", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec2->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec2->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec3 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec3", "HTCC Nphe vs beta for pip sector 3", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec3->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec3->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec4 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec4", "HTCC Nphe vs beta for pip sector 4", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec4->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec4->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec5 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec5", "HTCC Nphe vs beta for pip sector 5", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec5->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec5->GetYaxis()->SetTitle("#beta");
    hist_HTCC_Nphe_vs_beta_pip_sec6 = new TH2F("hist_HTCC_Nphe_vs_beta_pip_sec6", "HTCC Nphe vs beta for pip sector 6", 50, 0, 50, 700, 0.5, 1.2);
    hist_HTCC_Nphe_vs_beta_pip_sec6->GetXaxis()->SetTitle("Nphe");
    hist_HTCC_Nphe_vs_beta_pip_sec6->GetYaxis()->SetTitle("#beta");


    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    out->mkdir("Calorimeter_energy_depositions");
    out->cd("Calorimeter_energy_depositions");

    TH1F *hist_CAL_Edep_electron;
    TH1F *hist_CAL_Edep_proton;
    TH1F *hist_CAL_Edep_pip;
    TH1F *hist_CAL_Edep_pim;
    TH1F *hist_CAL_Edep_pip_diffsecele;
    TH1F *hist_CAL_Edep_pim_diffsecele;

    TH1F *hist_PCAL_Edep_electron;
    TH1F *hist_PCAL_Edep_proton;
    TH1F *hist_PCAL_Edep_pip;
    TH1F *hist_PCAL_Edep_pim;
    TH1F *hist_PCAL_Edep_pip_diffsecele;
    TH1F *hist_PCAL_Edep_pim_diffsecele;

    TH1F *hist_ECin_Edep_electron;
    TH1F *hist_ECin_Edep_proton;
    TH1F *hist_ECin_Edep_pip;
    TH1F *hist_ECin_Edep_pim;
    TH1F *hist_ECin_Edep_pip_diffsecele;
    TH1F *hist_ECin_Edep_pim_diffsecele;

    TH1F *hist_ECout_Edep_electron;
    TH1F *hist_ECout_Edep_proton;
    TH1F *hist_ECout_Edep_pip;
    TH1F *hist_ECout_Edep_pim;
    TH1F *hist_ECout_Edep_pip_diffsecele;
    TH1F *hist_ECout_Edep_pim_diffsecele;

    TH1F *hist_ECAL_Edep_electron;
    TH1F *hist_ECAL_Edep_proton;
    TH1F *hist_ECAL_Edep_pip;
    TH1F *hist_ECAL_Edep_pim;
    TH1F *hist_ECAL_Edep_pip_diffsecele;
    TH1F *hist_ECAL_Edep_pim_diffsecele;

    hist_CAL_Edep_electron = new TH1F("hist_CAL_Edep_electron", "Energy deposition of electrons in the complete Calorimeter", 500, 0, 2.5);
    hist_CAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_electron->GetYaxis()->SetTitle("counts");
    hist_CAL_Edep_proton = new TH1F("hist_CAL_Edep_proton", "Energy deposition of protons in the complete Calorimeter", 500, 0, 1);
    hist_CAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_proton->GetYaxis()->SetTitle("counts");
    hist_CAL_Edep_pip = new TH1F("hist_CAL_Edep_pip", "Energy deposition of pip in the complete Calorimeter", 600, 0, 1.5);
    hist_CAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_pip->GetYaxis()->SetTitle("counts");
    hist_CAL_Edep_pim = new TH1F("hist_CAL_Edep_pim", "Energy deposition of pim in the complete Calorimeter", 600, 0, 1.5);
    hist_CAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_pim->GetYaxis()->SetTitle("counts");
    hist_CAL_Edep_pip_diffsecele = new TH1F("hist_CAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in the complete Calorimeter", 600, 0, 1.5);
    hist_CAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
    hist_CAL_Edep_pim_diffsecele = new TH1F("hist_CAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in the complete Calorimeter", 600, 0, 1.5);
    hist_CAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_CAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

    hist_PCAL_Edep_electron = new TH1F("hist_PCAL_Edep_electron", "Energy deposition of electrons in PCAL", 400, 0, 2);
    hist_PCAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_electron->GetYaxis()->SetTitle("counts");
    hist_PCAL_Edep_proton = new TH1F("hist_PCAL_Edep_proton", "Energy deposition of protons in PCAL", 600, 0, 0.6);
    hist_PCAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_proton->GetYaxis()->SetTitle("counts");
    hist_PCAL_Edep_pip = new TH1F("hist_PCAL_Edep_pip", "Energy deposition of pip in PCAL", 500, 0, 1);
    hist_PCAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_pip->GetYaxis()->SetTitle("counts");
    hist_PCAL_Edep_pim = new TH1F("hist_PCAL_Edep_pim", "Energy deposition of pim in PCAL", 500, 0, 1);
    hist_PCAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_pim->GetYaxis()->SetTitle("counts");
    hist_PCAL_Edep_pip_diffsecele = new TH1F("hist_PCAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in PCAL", 500, 0, 1);
    hist_PCAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
    hist_PCAL_Edep_pim_diffsecele = new TH1F("hist_PCAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in PCAL", 500, 0, 1);
    hist_PCAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_PCAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

    hist_ECin_Edep_electron = new TH1F("hist_ECin_Edep_electron", "Energy deposition of electrons in ECin", 400, 0, 2);
    hist_ECin_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_electron->GetYaxis()->SetTitle("counts");
    hist_ECin_Edep_proton = new TH1F("hist_ECin_Edep_proton", "Energy deposition of protons in ECin", 500, 0, 1);
    hist_ECin_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_proton->GetYaxis()->SetTitle("counts");
    hist_ECin_Edep_pip = new TH1F("hist_ECin_Edep_pip", "Energy deposition of pip in ECin", 500, 0, 1);
    hist_ECin_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_pip->GetYaxis()->SetTitle("counts");
    hist_ECin_Edep_pim = new TH1F("hist_ECin_Edep_pim", "Energy deposition of pim in ECin", 500, 0, 1);
    hist_ECin_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_pim->GetYaxis()->SetTitle("counts");
    hist_ECin_Edep_pip_diffsecele = new TH1F("hist_ECin_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECin", 500, 0, 1);
    hist_ECin_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
    hist_ECin_Edep_pim_diffsecele = new TH1F("hist_ECin_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECin", 500, 0, 1);
    hist_ECin_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECin_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

    hist_ECout_Edep_electron = new TH1F("hist_ECout_Edep_electron", "Energy deposition of electrons in ECout", 400, 0, 2);
    hist_ECout_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_electron->GetYaxis()->SetTitle("counts");
    hist_ECout_Edep_proton = new TH1F("hist_ECout_Edep_proton", "Energy deposition of protons in ECout", 500, 0, 1);
    hist_ECout_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_proton->GetYaxis()->SetTitle("counts");
    hist_ECout_Edep_pip = new TH1F("hist_ECout_Edep_pip", "Energy deposition of pip in ECout", 500, 0, 1);
    hist_ECout_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_pip->GetYaxis()->SetTitle("counts");
    hist_ECout_Edep_pim = new TH1F("hist_ECout_Edep_pim", "Energy deposition of pim in ECout", 500, 0, 1);
    hist_ECout_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_pim->GetYaxis()->SetTitle("counts");
    hist_ECout_Edep_pip_diffsecele = new TH1F("hist_ECout_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECout", 500, 0, 1);
    hist_ECout_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
    hist_ECout_Edep_pim_diffsecele = new TH1F("hist_ECout_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECout", 500, 0, 1);
    hist_ECout_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECout_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");

    hist_ECAL_Edep_electron = new TH1F("hist_ECAL_Edep_electron", "Energy deposition of electrons in ECAL", 400, 0, 2);
    hist_ECAL_Edep_electron->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_electron->GetYaxis()->SetTitle("counts");
    hist_ECAL_Edep_proton = new TH1F("hist_ECAL_Edep_proton", "Energy deposition of protons in ECAL", 500, 0, 1);
    hist_ECAL_Edep_proton->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_proton->GetYaxis()->SetTitle("counts");
    hist_ECAL_Edep_pip = new TH1F("hist_ECAL_Edep_pip", "Energy deposition of pip in ECAL", 500, 0, 1);
    hist_ECAL_Edep_pip->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_pip->GetYaxis()->SetTitle("counts");
    hist_ECAL_Edep_pim = new TH1F("hist_ECAL_Edep_pim", "Energy deposition of pim in ECAL", 500, 0, 1);
    hist_ECAL_Edep_pim->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_pim->GetYaxis()->SetTitle("counts");
    hist_ECAL_Edep_pip_diffsecele = new TH1F("hist_ECAL_Edep_pip_diffsecele", "Energy deposition of pip (in a different sector than electrons) in ECAL", 500, 0, 1);
    hist_ECAL_Edep_pip_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_pip_diffsecele->GetYaxis()->SetTitle("counts");
    hist_ECAL_Edep_pim_diffsecele = new TH1F("hist_ECAL_Edep_pim_diffsecele", "Energy deposition of pim (in a different sector than electrons) in ECAL", 500, 0, 1);
    hist_ECAL_Edep_pim_diffsecele->GetXaxis()->SetTitle("E_{dep} /GeV");
    hist_ECAL_Edep_pim_diffsecele->GetYaxis()->SetTitle("counts");


    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    out->mkdir("statistics");
    out->cd("statistics");

    TH1F *hist_electron_count;
    TH1F *hist_proton_count;
    TH1F *hist_neutron_count;
    TH1F *hist_pip_count;
    TH1F *hist_pim_count;
    TH1F *hist_Kp_count;
    TH1F *hist_Km_count;
    TH1F *hist_photon_count;

    hist_electron_count = new TH1F("hist_electron_count", "electron count per event", 10, -0.5, 9.5);
    hist_electron_count->GetXaxis()->SetTitle("number of electrons per event");
    hist_electron_count->GetYaxis()->SetTitle("counts");
    hist_proton_count = new TH1F("hist_proton_count", "proton count per event", 10, -0.5, 9.5);
    hist_proton_count->GetXaxis()->SetTitle("number of protons per event");
    hist_proton_count->GetYaxis()->SetTitle("counts");
    hist_neutron_count = new TH1F("hist_neutron_count", "neutron count per event", 10, -0.5, 9.5);
    hist_neutron_count->GetXaxis()->SetTitle("number of neutrons per event");
    hist_neutron_count->GetYaxis()->SetTitle("counts");
    hist_pip_count = new TH1F("hist_pip_count", "pip count per event", 10, -0.5, 9.5);
    hist_pip_count->GetXaxis()->SetTitle("number of pips per event");
    hist_pip_count->GetYaxis()->SetTitle("counts");
    hist_pim_count = new TH1F("hist_pim_count", "pim count per event", 10, -0.5, 9.5);
    hist_pim_count->GetXaxis()->SetTitle("number of pims per event");
    hist_pim_count->GetYaxis()->SetTitle("counts");
    hist_Kp_count = new TH1F("hist_Kp_count", "Kp count per event", 10, -0.5, 9.5);
    hist_Kp_count->GetXaxis()->SetTitle("number of Kps per event");
    hist_Kp_count->GetYaxis()->SetTitle("counts");
    hist_Km_count = new TH1F("hist_Km_count", "Km count per event", 10, -0.5, 9.5);
    hist_Km_count->GetXaxis()->SetTitle("number of Kms per event");
    hist_Km_count->GetYaxis()->SetTitle("counts");
    hist_photon_count = new TH1F("hist_photon_count", "photon count per event", 10, -0.5, 9.5);
    hist_photon_count->GetXaxis()->SetTitle("number of photons per event");
    hist_photon_count->GetYaxis()->SetTitle("counts");

    out->mkdir("particles_charge");
    out->cd("particles_charge");

    TH1F *hist_particles_p;
    TH1F *hist_particles_theta;
    TH1F *hist_particles_phi;
    TH2F *hist_particles_p_vs_theta;
    TH2F *hist_particles_p_vs_phi;
    TH2F *hist_particles_theta_vs_phi;

    TH1F *hist_positives_p;
    TH1F *hist_positives_theta;
    TH1F *hist_positives_phi;
    TH2F *hist_positives_p_vs_theta;
    TH2F *hist_positives_p_vs_phi;
    TH2F *hist_positives_theta_vs_phi;

    TH1F *hist_negatives_p;
    TH1F *hist_negatives_theta;
    TH1F *hist_negatives_phi;
    TH2F *hist_negatives_p_vs_theta;
    TH2F *hist_negatives_p_vs_phi;
    TH2F *hist_negatives_theta_vs_phi;

    hist_particles_p = new TH1F("hist_particles_p", "particles momentum", 500, 0, Ebeam + 1);
    hist_particles_p->GetXaxis()->SetTitle("p /GeV");
    hist_particles_p->GetYaxis()->SetTitle("counts");
    hist_particles_theta = new TH1F("hist_particles_theta", "particles #Theta", 560, 0, 140);
    hist_particles_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_particles_theta->GetYaxis()->SetTitle("counts");
    hist_particles_phi = new TH1F("hist_particles_phi", "particles #phi", 180, -180, 180);
    hist_particles_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_particles_phi->GetYaxis()->SetTitle("counts");
    hist_particles_p_vs_theta = new TH2F("hist_particles_p_vs_theta", "particles p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_particles_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_particles_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_particles_p_vs_phi = new TH2F("hist_particles_p_vs_phi", "particles p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_particles_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_particles_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_particles_theta_vs_phi = new TH2F("hist_particles_theta_vs_phi", "particles #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_particles_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_particles_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_positives_p = new TH1F("hist_positives_p", "positives momentum", 500, 0, Ebeam + 1);
    hist_positives_p->GetXaxis()->SetTitle("p /GeV");
    hist_positives_p->GetYaxis()->SetTitle("counts");
    hist_positives_theta = new TH1F("hist_positives_theta", "positives #Theta", 560, 0, 140);
    hist_positives_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_positives_theta->GetYaxis()->SetTitle("counts");
    hist_positives_phi = new TH1F("hist_positives_phi", "positives #phi", 180, -180, 180);
    hist_positives_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_positives_phi->GetYaxis()->SetTitle("counts");
    hist_positives_p_vs_theta = new TH2F("hist_positives_p_vs_theta", "positives p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_positives_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_positives_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_positives_p_vs_phi = new TH2F("hist_positives_p_vs_phi", "positives p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_positives_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_positives_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_positives_theta_vs_phi = new TH2F("hist_positives_theta_vs_phi", "positives #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_positives_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_positives_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_negatives_p = new TH1F("hist_negatives_p", "negatives momentum", 500, 0, Ebeam + 1);
    hist_negatives_p->GetXaxis()->SetTitle("p /GeV");
    hist_negatives_p->GetYaxis()->SetTitle("counts");
    hist_negatives_theta = new TH1F("hist_negatives_theta", "negatives #Theta", 560, 0, 140);
    hist_negatives_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_negatives_theta->GetYaxis()->SetTitle("counts");
    hist_negatives_phi = new TH1F("hist_negatives_phi", "negatives #phi", 180, -180, 180);
    hist_negatives_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_negatives_phi->GetYaxis()->SetTitle("counts");
    hist_negatives_p_vs_theta = new TH2F("hist_negatives_p_vs_theta", "negatives p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_negatives_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_negatives_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_negatives_p_vs_phi = new TH2F("hist_negatives_p_vs_phi", "negatives p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_negatives_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_negatives_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_negatives_theta_vs_phi = new TH2F("hist_negatives_theta_vs_phi", "negatives #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_negatives_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_negatives_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    out->mkdir("particles_identified_histograms_all");
    out->cd("particles_identified_histograms_all");

    TH1F *hist_all_electron_p;
    TH1F *hist_all_electron_theta;
    TH1F *hist_all_electron_phi;
    TH2F *hist_all_electron_p_vs_theta;
    TH2F *hist_all_electron_p_vs_phi;
    TH2F *hist_all_electron_theta_vs_phi;

    TH1F *hist_all_proton_p;
    TH1F *hist_all_proton_theta;
    TH1F *hist_all_proton_phi;
    TH2F *hist_all_proton_p_vs_theta;
    TH2F *hist_all_proton_p_vs_phi;
    TH2F *hist_all_proton_theta_vs_phi;
    TH1F *hist_all_neutron_p;
    TH1F *hist_all_neutron_theta;
    TH1F *hist_all_neutron_phi;
    TH2F *hist_all_neutron_p_vs_theta;
    TH2F *hist_all_neutron_p_vs_phi;
    TH2F *hist_all_neutron_theta_vs_phi;

    TH1F *hist_all_pip_p;
    TH1F *hist_all_pip_theta;
    TH1F *hist_all_pip_phi;
    TH2F *hist_all_pip_p_vs_theta;
    TH2F *hist_all_pip_p_vs_phi;
    TH2F *hist_all_pip_theta_vs_phi;
    TH1F *hist_all_pim_p;
    TH1F *hist_all_pim_theta;
    TH1F *hist_all_pim_phi;
    TH2F *hist_all_pim_p_vs_theta;
    TH2F *hist_all_pim_p_vs_phi;
    TH2F *hist_all_pim_theta_vs_phi;

    TH1F *hist_all_Kp_p;
    TH1F *hist_all_Kp_theta;
    TH1F *hist_all_Kp_phi;
    TH2F *hist_all_Kp_p_vs_theta;
    TH2F *hist_all_Kp_p_vs_phi;
    TH2F *hist_all_Kp_theta_vs_phi;
    TH1F *hist_all_Km_p;
    TH1F *hist_all_Km_theta;
    TH1F *hist_all_Km_phi;
    TH2F *hist_all_Km_p_vs_theta;
    TH2F *hist_all_Km_p_vs_phi;
    TH2F *hist_all_Km_theta_vs_phi;

    TH1F *hist_all_photon_p;
    TH1F *hist_all_photon_theta;
    TH1F *hist_all_photon_phi;
    TH2F *hist_all_photon_p_vs_theta;
    TH2F *hist_all_photon_p_vs_phi;
    TH2F *hist_all_photon_theta_vs_phi;

    hist_all_electron_p = new TH1F("hist_all_electron_p", "electron momentum", 500, 0, Ebeam + 1);
    hist_all_electron_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_electron_p->GetYaxis()->SetTitle("counts");
    hist_all_electron_theta = new TH1F("hist_all_electron_theta", "electron #Theta", 200, 0, 50);
    hist_all_electron_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_electron_theta->GetYaxis()->SetTitle("counts");
    hist_all_electron_phi = new TH1F("hist_all_electron_phi", "electron #phi", 180, -180, 180);
    hist_all_electron_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_electron_phi->GetYaxis()->SetTitle("counts");
    hist_all_electron_p_vs_theta = new TH2F("hist_all_electron_p_vs_theta", "electron p vs #Theta", 200, 0, 50, 500, 0, Ebeam + 1);
    hist_all_electron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_electron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_electron_p_vs_phi = new TH2F("hist_all_electron_p_vs_phi", "electron p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_electron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_electron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_electron_theta_vs_phi = new TH2F("hist_all_electron_theta_vs_phi", "electron #Theta vs #phi", 180, -180, 180, 200, 0, 50);
    hist_all_electron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_electron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_proton_p = new TH1F("hist_all_proton_p", "proton momentum", 500, 0, Ebeam + 1);
    hist_all_proton_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_proton_p->GetYaxis()->SetTitle("counts");
    hist_all_proton_theta = new TH1F("hist_all_proton_theta", "proton #Theta", 560, 0, 140);
    hist_all_proton_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_proton_theta->GetYaxis()->SetTitle("counts");
    hist_all_proton_phi = new TH1F("hist_all_proton_phi", "proton #phi", 180, -180, 180);
    hist_all_proton_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_proton_phi->GetYaxis()->SetTitle("counts");
    hist_all_proton_p_vs_theta = new TH2F("hist_all_proton_p_vs_theta", "proton p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_proton_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_proton_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_proton_p_vs_phi = new TH2F("hist_all_proton_p_vs_phi", "proton p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_proton_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_proton_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_proton_theta_vs_phi = new TH2F("hist_all_proton_theta_vs_phi", "proton #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_proton_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_proton_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_neutron_p = new TH1F("hist_all_neutron_p", "neutron momentum", 500, 0, Ebeam + 1);
    hist_all_neutron_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_neutron_p->GetYaxis()->SetTitle("counts");
    hist_all_neutron_theta = new TH1F("hist_all_neutron_theta", "neutron #Theta", 560, 0, 140);
    hist_all_neutron_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_neutron_theta->GetYaxis()->SetTitle("counts");
    hist_all_neutron_phi = new TH1F("hist_all_neutron_phi", "neutron #phi", 180, -180, 180);
    hist_all_neutron_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_neutron_phi->GetYaxis()->SetTitle("counts");
    hist_all_neutron_p_vs_theta = new TH2F("hist_all_neutron_p_vs_theta", "neutron p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_neutron_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_neutron_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_neutron_p_vs_phi = new TH2F("hist_all_neutron_p_vs_phi", "neutron p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_neutron_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_neutron_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_neutron_theta_vs_phi = new TH2F("hist_all_neutron_theta_vs_phi", "neutron #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_neutron_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_neutron_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_pip_p = new TH1F("hist_all_pip_p", "#pi^{+} momentum", 500, 0, Ebeam + 1);
    hist_all_pip_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_pip_p->GetYaxis()->SetTitle("counts");
    hist_all_pip_theta = new TH1F("hist_all_pip_theta", "#pi^{+} #Theta", 560, 0, 140);
    hist_all_pip_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_pip_theta->GetYaxis()->SetTitle("counts");
    hist_all_pip_phi = new TH1F("hist_all_pip_phi", "#pi^{+} #phi", 180, -180, 180);
    hist_all_pip_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pip_phi->GetYaxis()->SetTitle("counts");
    hist_all_pip_p_vs_theta = new TH2F("hist_all_pip_p_vs_theta", "#pi^{+} p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_pip_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_pip_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_pip_p_vs_phi = new TH2F("hist_all_pip_p_vs_phi", "#pi^{+} p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_pip_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pip_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_pip_theta_vs_phi = new TH2F("hist_all_pip_theta_vs_phi", "#pi^{+} #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_pip_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pip_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_pim_p = new TH1F("hist_all_pim_p", "#pi^{-} momentum", 500, 0, Ebeam + 1);
    hist_all_pim_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_pim_p->GetYaxis()->SetTitle("counts");
    hist_all_pim_theta = new TH1F("hist_all_pim_theta", "#pi^{-} #Theta", 560, 0, 140);
    hist_all_pim_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_pim_theta->GetYaxis()->SetTitle("counts");
    hist_all_pim_phi = new TH1F("hist_all_pim_phi", "#pi^{-} #phi", 180, -180, 180);
    hist_all_pim_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pim_phi->GetYaxis()->SetTitle("counts");
    hist_all_pim_p_vs_theta = new TH2F("hist_all_pim_p_vs_theta", "#pi^{-} p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_pim_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_pim_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_pim_p_vs_phi = new TH2F("hist_all_pim_p_vs_phi", "#pi^{-} p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_pim_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pim_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_pim_theta_vs_phi = new TH2F("hist_all_pim_theta_vs_phi", "#pi^{-} #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_pim_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_pim_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_Kp_p = new TH1F("hist_all_Kp_p", "K^{+} momentum", 500, 0, Ebeam + 1);
    hist_all_Kp_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_Kp_p->GetYaxis()->SetTitle("counts");
    hist_all_Kp_theta = new TH1F("hist_all_Kp_theta", "K^{+} #Theta", 560, 0, 140);
    hist_all_Kp_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_Kp_theta->GetYaxis()->SetTitle("counts");
    hist_all_Kp_phi = new TH1F("hist_all_Kp_phi", "K^{+} #phi", 360, -180, 180);
    hist_all_Kp_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Kp_phi->GetYaxis()->SetTitle("counts");
    hist_all_Kp_p_vs_theta = new TH2F("hist_all_Kp_p_vs_theta", "K^{+} p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_Kp_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_Kp_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_Kp_p_vs_phi = new TH2F("hist_all_Kp_p_vs_phi", "K^{+} p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_Kp_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Kp_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_Kp_theta_vs_phi = new TH2F("hist_all_Kp_theta_vs_phi", "K^{+} #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_Kp_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Kp_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_Km_p = new TH1F("hist_all_Km_p", "K^{-} momentum", 500, 0, Ebeam + 1);
    hist_all_Km_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_Km_p->GetYaxis()->SetTitle("counts");
    hist_all_Km_theta = new TH1F("hist_all_Km_theta", "K^{-} #Theta", 560, 0, 140);
    hist_all_Km_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_Km_theta->GetYaxis()->SetTitle("counts");
    hist_all_Km_phi = new TH1F("hist_all_Km_phi", "K^{-} #phi", 180, -180, 180);
    hist_all_Km_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Km_phi->GetYaxis()->SetTitle("counts");
    hist_all_Km_p_vs_theta = new TH2F("hist_all_Km_p_vs_theta", "K^{-} p vs #Theta", 560, 0, 140, 500, 0, Ebeam + 1);
    hist_all_Km_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_Km_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_Km_p_vs_phi = new TH2F("hist_all_Km_p_vs_phi", "K^{-} p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_Km_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Km_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_Km_theta_vs_phi = new TH2F("hist_all_Km_theta_vs_phi", "K^{-} #Theta vs #phi", 180, -180, 180, 560, 0, 140);
    hist_all_Km_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_Km_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");

    hist_all_photon_p = new TH1F("hist_all_photon_p", "photon momentum", 500, 0, Ebeam + 1);
    hist_all_photon_p->GetXaxis()->SetTitle("p /GeV");
    hist_all_photon_p->GetYaxis()->SetTitle("counts");
    hist_all_photon_theta = new TH1F("hist_all_photon_theta", "photon #Theta", 200, 0, 50);
    hist_all_photon_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_photon_theta->GetYaxis()->SetTitle("counts");
    hist_all_photon_phi = new TH1F("hist_all_photon_phi", "photon #phi", 180, -180, 180);
    hist_all_photon_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_photon_phi->GetYaxis()->SetTitle("counts");
    hist_all_photon_p_vs_theta = new TH2F("hist_all_photon_p_vs_theta", "photon p vs #Theta", 200, 0, 50, 500, 0, Ebeam + 1);
    hist_all_photon_p_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_all_photon_p_vs_theta->GetYaxis()->SetTitle("p /GeV");
    hist_all_photon_p_vs_phi = new TH2F("hist_all_photon_p_vs_phi", "photon p vs #phi", 180, -180, 180, 500, 0, Ebeam + 1);
    hist_all_photon_p_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_photon_p_vs_phi->GetYaxis()->SetTitle("p /GeV");
    hist_all_photon_theta_vs_phi = new TH2F("hist_all_photon_theta_vs_phi", "photon #Theta vs #phi", 180, -180, 180, 200, 0, 50);
    hist_all_photon_theta_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_all_photon_theta_vs_phi->GetYaxis()->SetTitle("#Theta /deg");


    out->mkdir("vertex_plots");
    out->cd("vertex_plots");

    TH1F *hist_positive_vertex;
    TH1F *hist_negative_vertex;
    TH1F *hist_electron_vertex;
    TH1F *hist_proton_vertex;
    TH1F *hist_positive_vertex_sec[6];
    TH1F *hist_negative_vertex_sec[6];
    TH1F *hist_electron_vertex_sec[6];
    TH1F *hist_proton_vertex_sec[6];

    TH2F *hist_positive_vertex_vs_theta;
    TH2F *hist_negative_vertex_vs_theta;
    TH2F *hist_electron_vertex_vs_theta;
    TH2F *hist_proton_vertex_vs_theta;
    TH2F *hist_positive_vertex_vs_theta_sec[6];
    TH2F *hist_negative_vertex_vs_theta_sec[6];
    TH2F *hist_electron_vertex_vs_theta_sec[6];
    TH2F *hist_proton_vertex_vs_theta_sec[6];

    TH2F *hist_positive_vertex_vs_phi;
    TH2F *hist_negative_vertex_vs_phi;
    TH2F *hist_electron_vertex_vs_phi;
    TH2F *hist_proton_vertex_vs_phi;
    TH2F *hist_positive_vertex_vs_phi_sec[6];
    TH2F *hist_negative_vertex_vs_phi_sec[6];
    TH2F *hist_electron_vertex_vs_phi_sec[6];
    TH2F *hist_proton_vertex_vs_phi_sec[6];

    TH2F *hist_positive_vertex_vs_p;
    TH2F *hist_negative_vertex_vs_p;
    TH2F *hist_neutral_vertex_vs_p;
    TH2F *hist_electron_vertex_vs_p;
    TH2F *hist_proton_vertex_vs_p;
    TH2F *hist_positive_vertex_vs_p_sec[6];
    TH2F *hist_negative_vertex_vs_p_sec[6];
    TH2F *hist_neutral_vertex_vs_p_sec[6];
    TH2F *hist_electron_vertex_vs_p_sec[6];
    TH2F *hist_proton_vertex_vs_p_sec[6];

    hist_positive_vertex = new TH1F("hist_positive_vertex", "z vertex for particles with positive charge", 1000, -100, 100);
    hist_positive_vertex->GetXaxis()->SetTitle("v_{z} /cm");
    hist_positive_vertex->GetYaxis()->SetTitle("counts");
    hist_negative_vertex = new TH1F("hist_negative_vertex", "z vertex for particles with negative charge", 1000, -100, 100);
    hist_negative_vertex->GetXaxis()->SetTitle("v_{z} /cm");
    hist_negative_vertex->GetYaxis()->SetTitle("counts");
    hist_electron_vertex = new TH1F("hist_electron_vertex", "z vertex for electrons", 1000, -100, 100);
    hist_electron_vertex->GetXaxis()->SetTitle("v_{z} /cm");
    hist_electron_vertex->GetYaxis()->SetTitle("counts");
    hist_proton_vertex = new TH1F("hist_proton_vertex", "z vertex for protons", 1000, -100, 100);
    hist_proton_vertex->GetXaxis()->SetTitle("v_{z} /cm");
    hist_proton_vertex->GetYaxis()->SetTitle("counts");

    hist_positive_vertex_vs_theta = new TH2F("hist_positive_vertex_vs_theta", "z vertex vs #Theta for particles with positive charge", 140, 0, 140, 1000, -100, 100);
    hist_positive_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_positive_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
    hist_negative_vertex_vs_theta = new TH2F("hist_negative_vertex_vs_theta", "z vertex vs #Theta for particles with negative charge", 140, 0, 140, 1000, -100, 100);
    hist_negative_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_negative_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
    hist_electron_vertex_vs_theta = new TH2F("hist_electron_vertex_vs_theta", "z vertex vs #Theta for electrons", 140, 0, 140, 1000, -100, 100);
    hist_electron_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_electron_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");
    hist_proton_vertex_vs_theta = new TH2F("hist_proton_vertex_vs_theta", "z vertex vs #Theta for protons", 140, 0, 140, 1000, -100, 100);
    hist_proton_vertex_vs_theta->GetXaxis()->SetTitle("#Theta /deg");
    hist_proton_vertex_vs_theta->GetYaxis()->SetTitle("v_{z} /cm");

    hist_positive_vertex_vs_phi = new TH2F("hist_positive_vertex_vs_phi", "z vertex vs #phi for particles with positive charge", 140, 0, 140, 1000, -100, 100);
    hist_positive_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_positive_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
    hist_negative_vertex_vs_phi = new TH2F("hist_negative_vertex_vs_phi", "z vertex vs #phi for particles with negative charge", 140, 0, 140, 1000, -100, 100);
    hist_negative_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_negative_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
    hist_electron_vertex_vs_phi = new TH2F("hist_electron_vertex_vs_phi", "z vertex vs #phi for electrons", 140, 0, 140, 1000, -100, 100);
    hist_electron_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_electron_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");
    hist_proton_vertex_vs_phi = new TH2F("hist_proton_vertex_vs_phi", "z vertex vs #phi for protons", 140, 0, 140, 1000, -100, 100);
    hist_proton_vertex_vs_phi->GetXaxis()->SetTitle("#phi /deg");
    hist_proton_vertex_vs_phi->GetYaxis()->SetTitle("v_{z} /cm");

    hist_positive_vertex_vs_p = new TH2F("hist_positive_vertex_vs_p", "z vertex vs #p for particles with positive charge", 500, 0, Ebeam + 1, 1000, -100, 100);
    hist_positive_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV");
    hist_positive_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
    hist_negative_vertex_vs_p = new TH2F("hist_negative_vertex_vs_p", "z vertex vs #p for particles with negative charge", 500, 0, Ebeam + 1, 1000, -100, 100);
    hist_negative_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV");
    hist_negative_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
    hist_electron_vertex_vs_p = new TH2F("hist_electron_vertex_vs_p", "z vertex vs #p for electrons", 500, 0, Ebeam + 1, 1000, -100, 100);
    hist_electron_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV");
    hist_electron_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");
    hist_proton_vertex_vs_p = new TH2F("hist_proton_vertex_vs_p", "z vertex vs #p for protons", 500, 0, Ebeam + 1, 1000, -100, 100);
    hist_proton_vertex_vs_p->GetXaxis()->SetTitle("#p /GeV");
    hist_proton_vertex_vs_p->GetYaxis()->SetTitle("v_{z} /cm");

    for (Int_t i = 0; i < 6; i++)
    {

        sprintf(name, "hist_positive_vertex_sec%01d", i + 1);
        sprintf(title, "z vertex for particles with positive charge in sector %01d", i + 1);
        hist_positive_vertex_sec[i] = new TH1F(name, title, 1000, -100, 100);
        hist_positive_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm");
        hist_positive_vertex_sec[i]->GetYaxis()->SetTitle("counts");
        sprintf(name, "hist_negative_vertex_sec%01d", i + 1);
        sprintf(title, "z vertex for particles with negative charge in sector %01d", i + 1);
        hist_negative_vertex_sec[i] = new TH1F(name, title, 1000, -100, 100);
        hist_negative_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm");
        hist_negative_vertex_sec[i]->GetYaxis()->SetTitle("counts");
        sprintf(name, "hist_electron_vertex_sec%01d", i + 1);
        sprintf(title, "z vertex for electrons in sector %01d", i + 1);
        hist_electron_vertex_sec[i] = new TH1F(name, title, 1000, -100, 100);
        hist_electron_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm");
        hist_electron_vertex_sec[i]->GetYaxis()->SetTitle("counts");
        sprintf(name, "hist_proton_vertex_sec%01d", i + 1);
        sprintf(title, "z vertex for protons in sector %01d", i + 1);
        hist_proton_vertex_sec[i] = new TH1F(name, title, 1000, -100, 100);
        hist_proton_vertex_sec[i]->GetXaxis()->SetTitle("v_{z} /cm");
        hist_proton_vertex_sec[i]->GetYaxis()->SetTitle("counts");

        sprintf(name, "hist_positive_vertex_vs_theta_sec%01d", i + 1);
        sprintf(title, "z vertex vs #Theta for particles with positive charge in sector %01d", i + 1);
        hist_positive_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140, 1000, -100, 100);
        hist_positive_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg");
        hist_positive_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_negative_vertex_vs_theta_sec%01d", i + 1);
        sprintf(title, "z vertex vs #Theta for particles with negative charge in sector %01d", i + 1);
        hist_negative_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140, 1000, -100, 100);
        hist_negative_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg");
        hist_negative_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_electrons_vertex_vs_theta_sec%01d", i + 1);
        sprintf(title, "z vertex vs #Theta for electrons in sector %01d", i + 1);
        hist_electron_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140, 1000, -100, 100);
        hist_electron_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg");
        hist_electron_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_protons_vertex_vs_theta_sec%01d", i + 1);
        sprintf(title, "z vertex vs #Theta for protons in sector %01d", i + 1);
        hist_proton_vertex_vs_theta_sec[i] = new TH2F(name, title, 140, 0, 140, 1000, -100, 100);
        hist_proton_vertex_vs_theta_sec[i]->GetXaxis()->SetTitle("#Theta /deg");
        hist_proton_vertex_vs_theta_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");

        sprintf(name, "hist_positive_vertex_vs_phi_sec%01d", i + 1);
        sprintf(title, "z vertex vs #phi for particles with positive charge in sector %01d", i + 1);
        hist_positive_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180, 1000, -100, 100);
        hist_positive_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg");
        hist_positive_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_negative_vertex_vs_phi_sec%01d", i + 1);
        sprintf(title, "z vertex vs #phi for particles with negative charge in sector %01d", i + 1);
        hist_negative_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180, 1000, -100, 100);
        hist_negative_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg");
        hist_negative_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_electrons_vertex_vs_phi_sec%01d", i + 1);
        sprintf(title, "z vertex vs #phi for electrons in sector %01d", i + 1);
        hist_electron_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180, 1000, -100, 100);
        hist_electron_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg");
        hist_electron_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_protons_vertex_vs_phi_sec%01d", i + 1);
        sprintf(title, "z vertex vs #phi for protons in sector %01d", i + 1);
        hist_proton_vertex_vs_phi_sec[i] = new TH2F(name, title, 180, -180, 180, 1000, -100, 100);
        hist_proton_vertex_vs_phi_sec[i]->GetXaxis()->SetTitle("#phi /deg");
        hist_proton_vertex_vs_phi_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");

        sprintf(name, "hist_positive_vertex_vs_p_sec%01d", i + 1);
        sprintf(title, "z vertex vs p for particles with positive charge in sector %01d", i + 1);
        hist_positive_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam + 1, 1000, -100, 100);
        hist_positive_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV");
        hist_positive_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_negative_vertex_vs_p_sec%01d", i + 1);
        sprintf(title, "z vertex vs p for particles with negative charge in sector %01d", i + 1);
        hist_negative_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam + 1, 1000, -100, 100);
        hist_negative_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV");
        hist_negative_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_electrons_vertex_vs_p_sec%01d", i + 1);
        sprintf(title, "z vertex vs p for electrons in sector %01d", i + 1);
        hist_electron_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam + 1, 1000, -100, 100);
        hist_electron_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV");
        hist_electron_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
        sprintf(name, "hist_protons_vertex_vs_p_sec%01d", i + 1);
        sprintf(title, "z vertex vs p for protons in sector %01d", i + 1);
        hist_proton_vertex_vs_p_sec[i] = new TH2F(name, title, 500, 0, Ebeam + 1, 1000, -100, 100);
        hist_proton_vertex_vs_p_sec[i]->GetXaxis()->SetTitle("p /GeV");
        hist_proton_vertex_vs_p_sec[i]->GetYaxis()->SetTitle("v_{z} /cm");
    }

    out->mkdir("track_quality");
    out->cd("track_quality");

    TH1F *hist_FD_chi2;
    TH1F *hist_FD_chi2_NDF;
    TH1F *hist_FD_chi2_sec[6];
    TH1F *hist_FD_chi2_NDF_sec[6];
    TH1F *hist_CD_chi2;
    TH1F *hist_CD_chi2_NDF;

    hist_FD_chi2 = new TH1F("hist_FD_chi2", "chi2", 2000, 0, 1000);
    hist_FD_chi2->GetXaxis()->SetTitle("chi2");
    hist_FD_chi2->GetYaxis()->SetTitle("counts");
    hist_FD_chi2_NDF = new TH1F("hist_FD_chi2_NDF", "chi2/NDF", 1000, 0, 500);
    hist_FD_chi2_NDF->GetXaxis()->SetTitle("chi2/NDF");
    hist_FD_chi2_NDF->GetYaxis()->SetTitle("counts");

    for (Int_t i = 0; i < 6; i++)
    {
        sprintf(name, "hist_FD_chi2_sec%01d", i + 1);
        sprintf(title, "chi2 for sector %01d", i + 1);
        hist_FD_chi2_sec[i] = new TH1F(name, title, 2000, 0, 1000);
        hist_FD_chi2_sec[i]->GetXaxis()->SetTitle("chi2");
        hist_FD_chi2_sec[i]->GetYaxis()->SetTitle("counts");
        sprintf(name, "hist_FD_chi2_NDF_sec%01d", i + 1);
        sprintf(title, "chi2/NDF for sector %01d", i + 1);
        hist_FD_chi2_NDF_sec[i] = new TH1F(name, title, 1000, 0, 500);
        hist_FD_chi2_NDF_sec[i]->GetXaxis()->SetTitle("chi2/NDF");
        hist_FD_chi2_NDF_sec[i]->GetYaxis()->SetTitle("counts");
    }

    hist_CD_chi2 = new TH1F("hist_CD_chi2", "chi2", 2000, 0, 1000);
    hist_CD_chi2->GetXaxis()->SetTitle("chi2");
    hist_CD_chi2->GetYaxis()->SetTitle("counts");
    hist_CD_chi2_NDF = new TH1F("hist_CD_chi2_NDF", "chi2/NDF", 1000, 0, 500);
    hist_CD_chi2_NDF->GetXaxis()->SetTitle("chi2/NDF");
    hist_CD_chi2_NDF->GetYaxis()->SetTitle("counts");


    out->mkdir("additional_plots");
    out->cd("additional_plots");

    TH2F *hist_FTOF_phi_vs_sector_positives;
    TH2F *hist_FTOF_phi_vs_sector_negatives;

    hist_FTOF_phi_vs_sector_positives = new TH2F("FTOF_phi_vs_sector_positives", "#phi vs FTOF sector for particles with positive charge", 6, 0.5, 6.5, 180, -180, 180);
    hist_FTOF_phi_vs_sector_positives->GetXaxis()->SetTitle("FTOF sector");
    hist_FTOF_phi_vs_sector_positives->GetYaxis()->SetTitle("#phi /deg");

    hist_FTOF_phi_vs_sector_negatives = new TH2F("FTOF_phi_vs_sector_negatives", "#phi vs FTOF sector for particles with negative charge", 6, 0.5, 6.5, 180, -180, 180);
    hist_FTOF_phi_vs_sector_negatives->GetXaxis()->SetTitle("FTOF sector");
    hist_FTOF_phi_vs_sector_negatives->GetYaxis()->SetTitle("#phi /deg");

    out->mkdir("sampling_fraction");
    out->cd("sampling_fraction");

    TH2F *hist_sampling_fraction_vs_E_sec1;
    TH2F *hist_sampling_fraction_vs_E_sec2;
    TH2F *hist_sampling_fraction_vs_E_sec3;
    TH2F *hist_sampling_fraction_vs_E_sec4;
    TH2F *hist_sampling_fraction_vs_E_sec5;
    TH2F *hist_sampling_fraction_vs_E_sec6;

    TH2F *hist_sampling_fraction_vs_p_sec1;
    TH2F *hist_sampling_fraction_vs_p_sec2;
    TH2F *hist_sampling_fraction_vs_p_sec3;
    TH2F *hist_sampling_fraction_vs_p_sec4;
    TH2F *hist_sampling_fraction_vs_p_sec5;
    TH2F *hist_sampling_fraction_vs_p_sec6;

    hist_sampling_fraction_vs_E_sec1 = new TH2F("sampling_fraction_vs_E_sec1", "EC total sampling fraction versus E sec 1", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec1->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec1->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_E_sec2 = new TH2F("sampling_fraction_vs_E_sec2", "EC total sampling fraction versus E sec 2", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec2->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec2->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_E_sec3 = new TH2F("sampling_fraction_vs_E_sec3", "EC total sampling fraction versus E sec 3", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec3->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec3->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_E_sec4 = new TH2F("sampling_fraction_vs_E_sec4", "EC total sampling fraction versus E sec 4", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec4->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec4->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_E_sec5 = new TH2F("sampling_fraction_vs_E_sec5", "EC total sampling fraction versus E sec 5", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec5->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec5->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_E_sec6 = new TH2F("sampling_fraction_vs_E_sec6", "EC total sampling fraction versus E sec 6", 350, 0, 3.5, 300, 0, 0.6);
    hist_sampling_fraction_vs_E_sec6->GetXaxis()->SetTitle("E /GeV");
    hist_sampling_fraction_vs_E_sec6->GetYaxis()->SetTitle("E/p");

    hist_sampling_fraction_vs_p_sec1 = new TH2F("sampling_fraction_vs_p_sec1", "EC total sampling fraction versus p sec 1", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec1->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec1->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_p_sec2 = new TH2F("sampling_fraction_vs_p_sec2", "EC total sampling fraction versus p sec 2", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec2->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec2->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_p_sec3 = new TH2F("sampling_fraction_vs_p_sec3", "EC total sampling fraction versus p sec 3", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec3->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec3->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_p_sec4 = new TH2F("sampling_fraction_vs_p_sec4", "EC total sampling fraction versus p sec 4", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec4->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec4->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_p_sec5 = new TH2F("sampling_fraction_vs_p_sec5", "EC total sampling fraction versus p sec 5", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec5->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec5->GetYaxis()->SetTitle("E/p");
    hist_sampling_fraction_vs_p_sec6 = new TH2F("sampling_fraction_vs_p_sec6", "EC total sampling fraction versus p sec 6", 550, 0, 11, 300, 0, 0.6);
    hist_sampling_fraction_vs_p_sec6->GetXaxis()->SetTitle("p /GeV");
    hist_sampling_fraction_vs_p_sec6->GetYaxis()->SetTitle("E/p");
    

    out->mkdir("MC_LUND");
    out->cd("MC_LUND");

    TH1F *hist_LUND_pid;

    hist_LUND_pid = new TH1F("hist_LUND_pid", "generated PID", 6000, -3000, 3000);
    hist_LUND_pid->GetXaxis()->SetTitle("PID");
    hist_LUND_pid->GetYaxis()->SetTitle("counts");


    /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///  start of the event loop     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    neg_part_count = 0;
    pos_part_count = 0;
    neut_part_count = 0;

    cout << "Event Loop starting ... " << endl;
    cout << endl;

    if (use_qa_cuts)
    {
        c12.applyQA("latest");
        c12.db()->qadb_requireOkForAsymmetry(true);
    }

    double events = c12.getReader().getEntries();
    int count = 0;
    

    while (c12.next() && ((process_Events != -1 && count < process_Events) || (process_Events == -1 )))
    {

        count += 1;
        if(count % 10000 == 0){
          double percent = 100*count/events;
          cout << count << "  of  " << events << " events  (" << percent << ")" << endl;
        }

        std::vector<region_part_ptr> particles = c12.getDetParticles(); // particles is now a std::vector of particles for this event

        mcpar_ptr mcparticles = simulation ? c12.mcparts() : nullptr;
          
        
        bool qa_check = false;

        if (use_qa_cuts)
        {
            qa_check = true; //event is skipped by the c12 loop
        }
        else
        {
            qa_check = true;
        }

        if (qa_check == true)
        {

            /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// initalisatize the particle fourvectors and their components:
            ///

            TLorentzVector beam(0, 0, Ebeam, Ebeam);
            TLorentzVector target(0, 0, 0, 0.93827);

            NRUN = 0;
            NEVENT = 0;
            TYPE = 0;
            TRG = 0;
            Helic = 0;
            Helic_raw = 0;
            EVNTime = 0;
            BCG = 0;
            STTime = 0;
            RFTime = 0;

            helicity = 0;
            beam_charge = 0;
            eventnumber = c12.runconfig()->getEvent(); 

            p4_ele_px.clear();
            p4_prot_px.clear();
            p4_neutr_px.clear();
            p4_pip_px.clear();
            p4_pim_px.clear();
            p4_Kp_px.clear();
            p4_Km_px.clear();
            p4_phot_px.clear();
            p4_ele_py.clear();
            p4_prot_py.clear();
            p4_neutr_py.clear();
            p4_pip_py.clear();
            p4_pim_py.clear();
            p4_Kp_py.clear();
            p4_Km_py.clear();
            p4_phot_py.clear();
            p4_ele_pz.clear();
            p4_prot_pz.clear();
            p4_neutr_pz.clear();
            p4_pip_pz.clear();
            p4_pim_pz.clear();
            p4_Kp_pz.clear();
            p4_Km_pz.clear();
            p4_phot_pz.clear();
            p4_ele_E.clear();
            p4_prot_E.clear();
            p4_neutr_E.clear();
            p4_pip_E.clear();
            p4_pim_E.clear();
            p4_Kp_E.clear();
            p4_Km_E.clear();
            p4_phot_E.clear();
            p4_ele_chi2pid.clear();
            p4_prot_chi2pid.clear();
            p4_neutr_chi2pid.clear();
            p4_pip_chi2pid.clear();
            p4_pim_chi2pid.clear();
            p4_Kp_chi2pid.clear();
            p4_Km_chi2pid.clear();
            p4_ele_dcx1.clear(); p4_pip_dcx1.clear(); p4_pim_dcx1.clear(); 
            p4_ele_dcy1.clear(); p4_pip_dcy1.clear(); p4_pim_dcy1.clear(); 
            p4_ele_dcx2.clear(); p4_pip_dcx2.clear(); p4_pim_dcx2.clear();
            p4_ele_dcy2.clear(); p4_pip_dcy2.clear(); p4_pim_dcy2.clear(); 
            p4_ele_dcx3.clear(); p4_pip_dcx3.clear(); p4_pim_dcx3.clear(); 
            p4_ele_dcy3.clear(); p4_pip_dcy3.clear(); p4_pim_dcy3.clear(); 
            p4_prot_dcx1.clear(); p4_prot_dcy1.clear(); p4_prot_dcz1.clear(); 
            p4_ele_pcalv.clear(); p4_phot_pcalv.clear(); p4_ele_pcalw.clear(); p4_phot_pcalw.clear(); 
	    p4_ele_dcedge1.clear(); p4_ele_dcedge2.clear(); p4_ele_dcedge3.clear();
            p4_prot_dcedge1.clear(); p4_prot_dcedge2.clear(); p4_prot_dcedge3.clear();
            p4_pip_dcedge1.clear(); p4_pip_dcedge2.clear(); p4_pip_dcedge3.clear();
            p4_pim_dcedge1.clear(); p4_pim_dcedge2.clear(); p4_pim_dcedge3.clear();
            ele_det.clear();
            prot_det.clear();
            neutr_det.clear();
            pip_det.clear();
            pim_det.clear();
            Kp_det.clear();
            Km_det.clear();
            phot_det.clear();
            ele_sec.clear();
            prot_sec.clear();
            neutr_sec.clear();
            pip_sec.clear();
            pim_sec.clear();
            Kp_sec.clear();
            Km_sec.clear();
            pip_RICHbest.clear();
            pim_RICHbest.clear();
            Kp_RICHbest.clear();
            Km_RICHbest.clear();
            phot_sec.clear();
            p4_gen_ele_px.clear();
            p4_gen_prot_px.clear();
            p4_gen_neutr_px.clear();
            p4_gen_pip_px.clear();
            p4_gen_pim_px.clear();
            p4_gen_Kp_px.clear();
            p4_gen_Km_px.clear();
            p4_gen_phot_px.clear();
            p4_gen_ele_py.clear();
            p4_gen_prot_py.clear();
            p4_gen_neutr_py.clear();
            p4_gen_pip_py.clear();
            p4_gen_pim_py.clear();
            p4_gen_Kp_py.clear();
            p4_gen_Km_py.clear();
            p4_gen_phot_py.clear();
            p4_gen_ele_pz.clear();
            p4_gen_prot_pz.clear();
            p4_gen_neutr_pz.clear();
            p4_gen_pip_pz.clear();
            p4_gen_pim_pz.clear();
            p4_gen_Kp_pz.clear();
            p4_gen_Km_pz.clear();
            p4_gen_phot_pz.clear();
            p4_gen_ele_E.clear();
            p4_gen_prot_E.clear();
            p4_gen_neutr_E.clear();
            p4_gen_pip_E.clear();
            p4_gen_pim_E.clear();
            p4_gen_Kp_E.clear();
            p4_gen_Km_E.clear();
            p4_gen_phot_E.clear();
            p4_gen_pi0_px.clear();
            p4_gen_pi0_py.clear();
            p4_gen_pi0_pz.clear();
            p4_gen_pi0_E.clear();

            p4_Kp_ml.clear();
            p4_Kp_beta.clear();
            p4_Kp_LTCC_nphe.clear();
            p4_pip_ml.clear();
            p4_pip_beta.clear();
            p4_pip_LTCC_nphe.clear();
            p4_Km_ml.clear();
            p4_Km_beta.clear();
            p4_Km_LTCC_nphe.clear();
            p4_pim_ml.clear();
            p4_pim_beta.clear();
            p4_pim_LTCC_nphe.clear();
            
            p4_gen_ele_motherpid.clear(); 
            p4_gen_prot_motherpid.clear(); 
            p4_gen_neutr_motherpid.clear(); 
            p4_gen_pip_motherpid.clear(); 
            p4_gen_pim_motherpid.clear(); 
            p4_gen_Kp_motherpid.clear(); 
            p4_gen_Km_motherpid.clear(); 
            p4_gen_phot_motherpid.clear(); 
            p4_gen_pi0_motherpid.clear();

            e_count = 0;
            p_count = 0;
            n_count = 0;
            pip_count = 0;
            pim_count = 0;
            Kp_count = 0;
            Km_count = 0;
            g_count = 0;
            e_MCcount = 0;
            p_MCcount = 0;
            n_MCcount = 0;
            pip_MCcount = 0;
            pim_MCcount = 0;
            Kp_MCcount = 0;
            Km_MCcount = 0;
            g_MCcount = 0;

            MC_helicity = 0;
            MC_Npart = 0;
            MC_Ebeam = 0;
            MC_weight = 0;

            for (Int_t i = 0; i < BUFFER; i++)
            {

                p4_ele[i].SetPxPyPzE(0, 0, 0, 0);
                p4_ele_raw[i].SetPxPyPzE(0, 0, 0, 0);
                p4_prot[i].SetPxPyPzE(0, 0, 0, 0);
                p4_neutr[i].SetPxPyPzE(0, 0, 0, 0);
                p4_pip[i].SetPxPyPzE(0, 0, 0, 0);
                p4_pim[i].SetPxPyPzE(0, 0, 0, 0);
                p4_Kp[i].SetPxPyPzE(0, 0, 0, 0);
                p4_Km[i].SetPxPyPzE(0, 0, 0, 0);
                p4_phot[i].SetPxPyPzE(0, 0, 0, 0);

                e_ind[i] = -1;
                p_ind[i] = -1;
                n_ind[i] = -1;
                pip_ind[i] = -1;
                pim_ind[i] = -1;
                Kp_ind[i] = -1;
                Km_ind[i] = -1;
                g_ind[i] = -1;

                ele_detect[i] = 0;
                prot_detect[i] = 0;
                neutr_detect[i] = 0;
                pip_detect[i] = 0;
                pim_detect[i] = 0;
                Kp_detect[i] = 0;
                Km_detect[i] = 0;
                phot_detect[i] = 0;

                ele_sector[i] = 0;
                prot_sector[i] = 0;
                neutr_sector[i] = 0;
                pip_sector[i] = 0;
                pim_sector[i] = 0;
                Kp_sector[i] = 0;
                Km_sector[i] = 0;
                phot_sector[i] = 0;

                ele_chi2pid[i] = 0; 
                prot_chi2pid[i] = 0; 
                neutr_chi2pid[i] = 0; 
                pip_chi2pid[i] = 0; 
                pim_chi2pid[i] = 0; 
                Kp_chi2pid[i] = 0; 
                Km_chi2pid[i] = 0;
                
                ele_dcx1[i] = 0; pip_dcx1[i] = 0; pim_dcx1[i] = 0; 
                ele_dcy1[i] = 0; pip_dcy1[i] = 0; pim_dcy1[i] = 0; 
                ele_dcx2[i] = 0; pip_dcx2[i] = 0; pim_dcx2[i] = 0;
                ele_dcy2[i] = 0; pip_dcy2[i] = 0; pim_dcy2[i] = 0;
                ele_dcx3[i] = 0; pip_dcx3[i] = 0; pim_dcx3[i] = 0;
                ele_dcy3[i] = 0; pip_dcy3[i] = 0; pim_dcy3[i] = 0;
                prot_dcx1[i] = 0; prot_dcy1[i] = 0; prot_dcz1[i] = 0; 
                ele_pcalv[i] = 0; phot_pcalv[i] = 0; ele_pcalw[i] = 0; phot_pcalw[i] = 0;
                ele_dcedge1[i] = 0; prot_dcedge1[i] = 0; pip_dcedge1[i] = 0; pim_dcedge1[i] = 0;
                ele_dcedge2[i] = 0; prot_dcedge1[i] = 0; pip_dcedge2[i] = 0; pim_dcedge2[i] = 0;
                ele_dcedge3[i] = 0; prot_dcedge1[i] = 0; pip_dcedge3[i] = 0; pim_dcedge3[i] = 0;

                e_vx[i] = 0;
                e_vy[i] = 0;
                e_vz[i] = 0;
                e_beta[i] = 0, e_FTOF_sec[i] = -1, e_PCAL_sec[i] = -1;
                p_vx[i] = 0;
                p_vy[i] = 0;
                p_vz[i] = 0;
                p_beta[i] = 0, p_FTOF_sec[i] = -1, p_PCAL_sec[i] = -1;
                n_vx[i] = 0;
                n_vy[i] = 0;
                n_vz[i] = 0;
                n_beta[i] = 0;
                pip_vx[i] = 0;
                pip_vy[i] = 0;
                pip_vz[i] = 0;
                pip_beta[i] = 0, pip_FTOF_sec[i] = -1;
                pim_vx[i] = 0;
                pim_vy[i] = 0;
                pim_vz[i] = 0;
                pim_beta[i] = 0, pim_FTOF_sec[i] = -1;
                Kp_vx[i] = 0;
                Kp_vy[i] = 0;
                Kp_vz[i] = 0;
                Kp_beta[i] = 0, Kp_FTOF_sec[i] = -1;
                Km_vx[i] = 0;
                Km_vy[i] = 0;
                Km_vz[i] = 0;
                Km_beta[i] = 0, Km_FTOF_sec[i] = -1;
                g_vx[i] = 0;
                g_vy[i] = 0;
                g_vz[i] = 0;
                g_sec[i] = 0;

                electron_sector_cut[i] = false;
            }

            for (Int_t i = 0; i < BUFFER; i++)
            {

                part_px[i] = 0;
                part_py[i] = 0;
                part_pz[i] = 0;
                part_p[i] = 0;
                part_beta[i] = 0;
                part_vx[i] = 0;
                part_vy[i] = 0;
                part_vz[i] = 0;
                partMC_mother[i] = 0;
                part_charge[i] = 0;
                part_pid[i] = 0;
                part_theta[i] = 0;
                part_phi[i] = 0;
                part_status[i] = 0;
                part_chi2pid[i] = 0;

                partMC_px[i] = 0;
                partMC_py[i] = 0;
                partMC_pz[i] = 0;
                partMC_p[i] = 0;
                partMC_theta[i] = 0;
                partMC_phi[i] = 0;
                partMC_vx[i] = 0;
                partMC_vy[i] = 0;
                partMC_vz[i] = 0;
                partMC_pid[i] = 0;
                partLUND_mass[i] = 0;
                partLUND_E[i] = 0;
                partLUND_px[i] = 0;
                partLUND_py[i] = 0;
                partLUND_pz[i] = 0;
                partLUND_p[i] = 0;
                partLUND_theta[i] = 0;
                partLUND_phi[i] = 0;
                partLUND_vx[i] = 0;
                partLUND_vy[i] = 0;
                partLUND_vz[i] = 0;
                partLUND_pid[i] = 0;
                partLUND_mother[i] = 0;

                part_Cal_PCAL_sector[i] = 0;
                part_Cal_ECin_sector[i] = 0;
                part_Cal_ECout_sector[i] = 0;
                part_Cal_PCAL_energy[i] = 0;
                part_Cal_ECin_energy[i] = 0;
                part_Cal_ECout_energy[i] = 0;
                part_Cal_energy_total[i] = 0;
                part_Cal_PCAL_time[i] = 0;
                part_Cal_ECin_time[i] = 0;
                part_Cal_ECout_time[i] = 0;
                part_Cal_PCAL_path[i] = 0;
                part_Cal_ECin_path[i] = 0;
                part_Cal_ECout_path[i] = 0;
                part_Cal_PCAL_x[i] = 0;
                part_Cal_PCAL_y[i] = 0;
                part_Cal_PCAL_z[i] = 0;
                part_Cal_ECin_x[i] = 0;
                part_Cal_ECin_y[i] = 0;
                part_Cal_ECin_z[i] = 0;
                part_Cal_ECout_x[i] = 0;
                part_Cal_ECout_y[i] = 0;
                part_Cal_ECout_z[i] = 0;
                part_Cal_PCAL_lu[i] = 0;
                part_Cal_PCAL_lv[i] = 0;
                part_Cal_PCAL_lw[i] = 0;
                part_Cal_ECin_lu[i] = 0;
                part_Cal_ECin_lv[i] = 0;
                part_Cal_ECin_lw[i] = 0;
                part_Cal_ECout_lu[i] = 0;
                part_Cal_ECout_lv[i] = 0;
                part_Cal_ECout_lw[i] = 0;

                part_CC_HTCC_sector[i] = 0;
                part_CC_HTCC_nphe[i] = 0;
                part_CC_HTCC_time[i] = 0;
                part_CC_HTCC_path[i] = 0;
                part_CC_HTCC_theta[i] = 0;
                part_CC_HTCC_phi[i] = 0;

                part_CC_LTCC_sector[i] = 0;
                part_CC_LTCC_nphe[i] = 0;
                part_CC_LTCC_time[i] = 0;
                part_CC_LTCC_path[i] = 0;
                part_CC_LTCC_theta[i] = 0;
                part_CC_LTCC_phi[i] = 0;

                part_FTOF_sector_layer1[i] = 0;
                part_FTOF_sector_layer2[i] = 0;
                part_FTOF_sector_layer3[i] = 0;
                part_FTOF_component_layer1[i] = 0;
                part_FTOF_component_layer2[i] = 0;
                part_FTOF_component_layer3[i] = 0;
                part_FTOF_energy[i] = 0;
                part_FTOF_time[i] = 0;
                part_FTOF_path[i] = 0;
                part_FTOF_energy_layer1[i] = 0;
                part_FTOF_time_layer1[i] = 0;
                part_FTOF_path_layer1[i] = 0;
                part_FTOF_energy_layer3[i] = 0;
                part_FTOF_time_layer3[i] = 0;
                part_FTOF_path_layer3[i] = 0;
                part_FTOF_layer[i] = 0;

                part_CTOF_component[i] = 0;
                part_CTOF_energy[i] = 0;
                part_CTOF_time[i] = 0;
                part_CTOF_path[i] = 0;

                part_CND_component[i] = 0;
                part_CND_energy[i] = 0;
                part_CND_time[i] = 0;
                part_CND_path[i] = 0;

                part_FT_energy[i] = 0;
                part_FT_radius[i] = 0;
                part_FTHODO_time[i] = 0;
                part_FTHODO_path[i] = 0;
                part_FTCAL_time[i] = 0;
                part_FTCAL_path[i] = 0;
                part_FTCAL_x[i] = 0;
                part_FTCAL_y[i] = 0;
                part_FTCAL_z[i] = 0;
                part_FTTRK_x[i] = 0;
                part_FTTRK_y[i] = 0;
                part_FTTRK_z[i] = 0;
                part_FTHODO_x[i] = 0;
                part_FTHODO_y[i] = 0;
                part_FTHODO_z[i] = 0;

                part_DC_index[i] = 0;
                part_DC_sector[i] = 0;
                part_DC_Track_chi2[i] = 0;
                part_DC_Track_NDF[i] = 0;
                part_DC_Track_status[i] = 0;
                part_DC_c1x[i] = 0;
                part_DC_c1y[i] = 0;
                part_DC_c1z[i] = 0;
                part_DC_c2x[i] = 0;
                part_DC_c2y[i] = 0;
                part_DC_c2z[i] = 0;
                part_DC_c3x[i] = 0;
                part_DC_c3y[i] = 0;
                part_DC_c3z[i] = 0;

                part_DC_edge1[i] = 0; 
                part_DC_edge3[i] = 0; 
                part_DC_edge3[i] = 0;

                part_RICH_best[i] = 0;

                FD_eid_default_PID_check[i] = false;
                FD_eid_charge_check[i] = false;
                FD_eid_EC_outer_vs_EC_inner_check[i] = false;
                FD_eid_EC_sampling_fraction_check[i] = false;
                FD_eid_EC_hit_position_fiducial_check[i] = false;
                FD_eid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_eid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_eid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_eid_DC_z_vertex_check[i] = false;
                FD_eid_all_check[i] = false;

                FD_protid_default_PID_check[i] = false;
                FD_protid_charge_check[i] = false;
                FD_protid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_protid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_protid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_protid_beta_check[i] = false;
                FD_protid_delta_beta_check[i] = false;
                FD_protid_tofmass_check[i] = false;
                FD_protid_maximum_probability_check[i] = false;
                FD_protid_delta_vz_check[i] = false;
                FD_protid_all_check[i] = false;
                FD_neutrid_default_PID_check[i] = false;
                FD_neutrid_charge_check[i] = false;
                FD_neutrid_beta_check[i] = false;
                FD_neutrid_delta_beta_check[i] = false;
                FD_neutrid_tofmass_check[i] = false;
                FD_neutrid_delta_vz_check[i] = false;
                FD_neutrid_all_check[i] = false;
                FD_pipid_default_PID_check[i] = false;
                FD_pipid_charge_check[i] = false;
                FD_pipid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_pipid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_pipid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_pipid_beta_check[i] = false;
                FD_pipid_delta_beta_check[i] = false;
                FD_pipid_tofmass_check[i] = false;
                FD_pipid_maximum_probability_check[i] = false;
                FD_pipid_delta_vz_check[i] = false;
                FD_pipid_all_check[i] = false;
                FD_pimid_default_PID_check[i] = false;
                FD_pimid_charge_check[i] = false;
                FD_pimid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_pimid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_pimid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_pimid_beta_check[i] = false;
                FD_pimid_delta_beta_check[i] = false;
                FD_pimid_tofmass_check[i] = false;
                FD_pimid_maximum_probability_check[i] = false;
                FD_pimid_delta_vz_check[i] = false;
                FD_pimid_all_check[i] = false;
                FD_Kpid_default_PID_check[i] = false;
                FD_Kpid_charge_check[i] = false;
                FD_Kpid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_Kpid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_Kpid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_Kpid_beta_check[i] = false;
                FD_Kpid_delta_beta_check[i] = false;
                FD_Kpid_tofmass_check[i] = false;
                FD_Kpid_maximum_probability_check[i] = false;
                FD_Kpid_delta_vz_check[i] = false;
                FD_Kpid_all_check[i] = false;
                FD_Kmid_default_PID_check[i] = false;
                FD_Kmid_charge_check[i] = false;
                FD_Kmid_DC_hit_position_region1_fiducial_check[i] = false;
                FD_Kmid_DC_hit_position_region2_fiducial_check[i] = false;
                FD_Kmid_DC_hit_position_region3_fiducial_check[i] = false;
                FD_Kmid_beta_check[i] = false;
                FD_Kmid_delta_beta_check[i] = false;
                FD_Kmid_tofmass_check[i] = false;
                FD_Kmid_maximum_probability_check[i] = false;
                FD_Kmid_delta_vz_check[i] = false;
                FD_Kmid_all_check[i] = false;

                FD_photid_default_PID_check[i] = false;
                FD_photid_charge_check[i] = false;
                FD_photid_beta_check[i] = false;
                FD_photid_EC_sampling_fraction_check[i] = false;
                FD_photid_EC_hit_position_fiducial_check[i] = false;
                FD_photid_all_check[i] = false;

                // FT

                FT_eid_charge_check[i] = false;
                FT_eid_PID_check[i] = false;
                FT_eid_FTCAL_fiducial_check[i] = false;
                FT_eid_FTTRK_fiducial_check[i] = false;
                FT_eid_FTHODO_fiducial_check[i] = false;
                FT_eid_energy_vs_radius_check[i] = false;
                FT_eid_all_check[i] = false;

                FT_photid_charge_check[i] = false;
                FT_photid_PID_check[i] = false;
                FT_photid_FTCAL_fiducial_check[i] = false;
                FT_photid_beta_check[i] = false;
                FT_photid_all_check[i] = false;

                // CD

                CD_protid_default_PID_check[i] = false;
                CD_protid_charge_check[i] = false;
                CD_protid_beta_check[i] = false;
                CD_protid_maximum_probability_check[i] = false;
                CD_protid_delta_vz_check[i] = false;
                CD_protid_all_check[i] = false;

                CD_neutrid_default_PID_check[i] = false;
                CD_neutrid_charge_check[i] = false;
                CD_neutrid_beta_check[i] = false;
                CD_neutrid_maximum_probability_check[i] = false;
                CD_neutrid_delta_vz_check[i] = false;
                CD_neutrid_all_check[i] = false;

                CD_pipid_default_PID_check[i] = false;
                CD_pipid_charge_check[i] = false;
                CD_pipid_beta_check[i] = false;
                CD_pipid_maximum_probability_check[i] = false;
                CD_pipid_delta_vz_check[i] = false;
                CD_pipid_all_check[i] = false;

                CD_pimid_default_PID_check[i] = false;
                CD_pimid_charge_check[i] = false;
                CD_pimid_beta_check[i] = false;
                CD_pimid_maximum_probability_check[i] = false;
                CD_pimid_delta_vz_check[i] = false;
                CD_pimid_all_check[i] = false;

                CD_Kpid_default_PID_check[i] = false;
                CD_Kpid_charge_check[i] = false;
                CD_Kpid_beta_check[i] = false;
                CD_Kpid_maximum_probability_check[i] = false;
                CD_Kpid_delta_vz_check[i] = false;
                CD_Kpid_all_check[i] = false;

                CD_Kmid_default_PID_check[i] = false;
                CD_Kmid_charge_check[i] = false;
                CD_Kmid_beta_check[i] = false;
                CD_Kmid_maximum_probability_check[i] = false;
                CD_Kmid_delta_vz_check[i] = false;
                CD_Kmid_all_check[i] = false;
            }


            /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //  get event properties and assign particles to variables
            ///

            get_event_properties(c12.event(), simulation? c12.mcevent() : nullptr);

            helicity = Helic;

            assign_particles(particles, mcparticles);

            if (simulation == true && MC_helicity != 0)
            {
                helicity = MC_helicity;
            }

            /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //  do particle ID an momentum correction
            ///

            if (simulation == true)
            {
                write_gen(run);
            }

            select_electron(run); // first select good electrons

            if (e_count > 0)
            { // basic trigger condition: at least one good electron

                // continue to select other particles
                select_proton(run);
                select_neutron(run);
                select_pip(run);
                select_pim(run);
                select_Kplus(run);
                select_Kminus(run);
                select_photon(run);


                // fill output tree variables:

                fill_output_vector_electron();
                fill_output_vector_proton();
                fill_output_vector_neutron();
                fill_output_vector_pip();
                fill_output_vector_pim();
                fill_output_vector_Kplus();
                fill_output_vector_Kminus();
                fill_output_vector_photon();

            } // end of basic trigger condition

            // write particles to the tree:
            

            out_tree.Fill();


            /// ////////////////////////////////////////////////////////////////////////////////////////
            /// fill the histograms
            /// ////////////////////////////////////////////////////////////////////////////////////////

            if (e_count > 0)
            {

                hist_beam_charge->Fill(beam_charge);
                hist_helicity->Fill(helicity);
                hist_event_starttime->Fill(STTime);
                hist_RF_time->Fill(RFTime);


                // fill particle distributions:

                for (Int_t i = 0; i < BUFFER; i++)
                {

                    if (part_status[i] > 0)
                        hist_status->Fill(part_status[i]);

                    if (part_p[i] > 0)
                        hist_particles_p->Fill(part_p[i]);
                    if (part_theta[i] > 0)
                        hist_particles_theta->Fill(part_theta[i] * 180 / Pival);
                    if (part_phi[i] != 0)
                        hist_particles_phi->Fill(part_phi[i] * 180 / Pival);
                    if (part_theta[i] > 0 && part_p[i] > 0)
                        hist_particles_p_vs_theta->Fill(part_theta[i] * 180 / Pival, part_p[i]);
                    if (part_phi[i] != 0 && part_p[i] > 0)
                        hist_particles_p_vs_phi->Fill(part_phi[i] * 180 / Pival, part_p[i]);
                    if (part_theta[i] > 0 && part_phi[i] != 0)
                        hist_particles_theta_vs_phi->Fill(part_phi[i] * 180 / Pival, part_theta[i] * 180 / Pival);

                    if (part_charge[i] == +1)
                    {
                        if (part_p[i] > 0)
                            hist_positives_p->Fill(part_p[i]);
                        if (part_theta[i] > 0)
                            hist_positives_theta->Fill(part_theta[i] * 180 / Pival);
                        if (part_phi[i] != 0)
                            hist_positives_phi->Fill(part_phi[i] * 180 / Pival);
                        if (part_theta[i] > 0 && part_p[i] > 0)
                            hist_positives_p_vs_theta->Fill(part_theta[i] * 180 / Pival, part_p[i]);
                        if (part_phi[i] != 0 && part_p[i] > 0)
                            hist_positives_p_vs_phi->Fill(part_phi[i] * 180 / Pival, part_p[i]);
                        if (part_theta[i] > 0 && part_phi[i] != 0)
                            hist_positives_theta_vs_phi->Fill(part_phi[i] * 180 / Pival, part_theta[i] * 180 / Pival);
                    }
                    if (part_charge[i] == -1)
                    {
                        if (part_p[i] > 0)
                            hist_negatives_p->Fill(part_p[i]);
                        if (part_theta[i] > 0)
                            hist_negatives_theta->Fill(part_theta[i] * 180 / Pival);
                        if (part_phi[i] != 0)
                            hist_negatives_phi->Fill(part_phi[i] * 180 / Pival);
                        if (part_theta[i] > 0 && part_p[i] > 0)
                            hist_negatives_p_vs_theta->Fill(part_theta[i] * 180 / Pival, part_p[i]);
                        if (part_phi[i] != 0 && part_p[i] > 0)
                            hist_negatives_p_vs_phi->Fill(part_phi[i] * 180 / Pival, part_p[i]);
                        if (part_theta[i] > 0 && part_phi[i] != 0)
                            hist_negatives_theta_vs_phi->Fill(part_phi[i] * 180 / Pival, part_theta[i] * 180 / Pival);
                    }
                }


                for (Int_t i = 0; i < BUFFER; i++)
                {
                    if (p4_ele[i].P() > 0)
                        hist_all_electron_p->Fill(p4_ele[i].P());
                    if (p4_ele[i].Theta() > 0)
                        hist_all_electron_theta->Fill(p4_ele[i].Theta() * 180 / Pival);
                    if (p4_ele[i].Phi() != 0)
                        hist_all_electron_phi->Fill(p4_ele[i].Phi() * 180 / Pival);
                    if (p4_ele[i].P() > 0)
                        hist_all_electron_p_vs_theta->Fill(p4_ele[i].Theta() * 180 / Pival, p4_ele[i].P());
                    if (p4_ele[i].P() > 0)
                        hist_all_electron_p_vs_phi->Fill(p4_ele[i].Phi() * 180 / Pival, p4_ele[i].P());
                    if (p4_ele[i].Theta() > 0 && p4_ele[i].Phi() != 0)
                        hist_all_electron_theta_vs_phi->Fill(p4_ele[i].Phi() * 180 / Pival, p4_ele[i].Theta() * 180 / Pival);

                    if (p4_prot[i].P() > 0)
                        hist_all_proton_p->Fill(p4_prot[i].P());
                    if (p4_prot[i].Theta() > 0)
                        hist_all_proton_theta->Fill(p4_prot[i].Theta() * 180 / Pival);
                    if (p4_prot[i].Phi() != 0)
                        hist_all_proton_phi->Fill(p4_prot[i].Phi() * 180 / Pival);
                    if (p4_prot[i].P() > 0)
                        hist_all_proton_p_vs_theta->Fill(p4_prot[i].Theta() * 180 / Pival, p4_prot[i].P());
                    if (p4_prot[i].P() > 0)
                        hist_all_proton_p_vs_phi->Fill(p4_prot[i].Phi() * 180 / Pival, p4_prot[i].P());
                    if (p4_prot[i].Theta() > 0 && p4_prot[i].Phi() != 0)
                        hist_all_proton_theta_vs_phi->Fill(p4_prot[i].Phi() * 180 / Pival, p4_prot[i].Theta() * 180 / Pival);

                    if (p4_neutr[i].P() > 0)
                        hist_all_neutron_p->Fill(p4_neutr[i].P());
                    if (p4_neutr[i].Theta() > 0)
                        hist_all_neutron_theta->Fill(p4_neutr[i].Theta() * 180 / Pival);
                    if (p4_neutr[i].Phi() != 0)
                        hist_all_neutron_phi->Fill(p4_neutr[i].Phi() * 180 / Pival);
                    if (p4_neutr[i].P() > 0)
                        hist_all_neutron_p_vs_theta->Fill(p4_neutr[i].Theta() * 180 / Pival, p4_neutr[i].P());
                    if (p4_neutr[i].P() > 0)
                        hist_all_neutron_p_vs_phi->Fill(p4_neutr[i].Phi() * 180 / Pival, p4_neutr[i].P());
                    if (p4_neutr[i].Theta() > 0 && p4_neutr[i].Phi() != 0)
                        hist_all_neutron_theta_vs_phi->Fill(p4_neutr[i].Phi() * 180 / Pival, p4_neutr[i].Theta() * 180 / Pival);

                    if (p4_pip[i].P() > 0)
                        hist_all_pip_p->Fill(p4_pip[i].P());
                    if (p4_pip[i].Theta() > 0)
                        hist_all_pip_theta->Fill(p4_pip[i].Theta() * 180 / Pival);
                    if (p4_pip[i].Phi() != 0)
                        hist_all_pip_phi->Fill(p4_pip[i].Phi() * 180 / Pival);
                    if (p4_pip[i].P() > 0)
                        hist_all_pip_p_vs_theta->Fill(p4_pip[i].Theta() * 180 / Pival, p4_pip[i].P());
                    if (p4_pip[i].P() > 0)
                        hist_all_pip_p_vs_phi->Fill(p4_pip[i].Phi() * 180 / Pival, p4_pip[i].P());
                    if (p4_pip[i].Theta() > 0 && p4_pip[i].Phi() != 0)
                        hist_all_pip_theta_vs_phi->Fill(p4_pip[i].Phi() * 180 / Pival, p4_pip[i].Theta() * 180 / Pival);

                    if (p4_pim[i].P() > 0)
                        hist_all_pim_p->Fill(p4_pim[i].P());
                    if (p4_pim[i].Theta() > 0)
                        hist_all_pim_theta->Fill(p4_pim[i].Theta() * 180 / Pival);
                    if (p4_pim[i].Phi() != 0)
                        hist_all_pim_phi->Fill(p4_pim[i].Phi() * 180 / Pival);
                    if (p4_pim[i].P() > 0)
                        hist_all_pim_p_vs_theta->Fill(p4_pim[i].Theta() * 180 / Pival, p4_pim[i].P());
                    if (p4_pim[i].P() > 0)
                        hist_all_pim_p_vs_phi->Fill(p4_pim[i].Phi() * 180 / Pival, p4_pim[i].P());
                    if (p4_pim[i].Theta() > 0 && p4_pim[i].Phi() != 0)
                        hist_all_pim_theta_vs_phi->Fill(p4_pim[i].Phi() * 180 / Pival, p4_pim[i].Theta() * 180 / Pival);

                    if (p4_Kp[i].P() > 0)
                        hist_all_Kp_p->Fill(p4_Kp[i].P());
                    if (p4_Kp[i].Theta() > 0)
                        hist_all_Kp_theta->Fill(p4_Kp[i].Theta() * 180 / Pival);
                    if (p4_Kp[i].Phi() != 0)
                        hist_all_Kp_phi->Fill(p4_Kp[i].Phi() * 180 / Pival);
                    if (p4_Kp[i].P() > 0)
                        hist_all_Kp_p_vs_theta->Fill(p4_Kp[i].Theta() * 180 / Pival, p4_Kp[i].P());
                    if (p4_Kp[i].P() > 0)
                        hist_all_Kp_p_vs_phi->Fill(p4_Kp[i].Phi() * 180 / Pival, p4_Kp[i].P());
                    if (p4_Kp[i].Theta() > 0 && p4_Kp[i].Phi() != 0)
                        hist_all_Kp_theta_vs_phi->Fill(p4_Kp[i].Phi() * 180 / Pival, p4_Kp[i].Theta() * 180 / Pival);

                    if (p4_Km[i].P() > 0)
                        hist_all_Km_p->Fill(p4_Km[i].P());
                    if (p4_Km[i].Theta() > 0)
                        hist_all_Km_theta->Fill(p4_Km[i].Theta() * 180 / Pival);
                    if (p4_Km[i].Phi() != 0)
                        hist_all_Km_phi->Fill(p4_Km[i].Phi() * 180 / Pival);
                    if (p4_Km[i].P() > 0)
                        hist_all_Km_p_vs_theta->Fill(p4_Km[i].Theta() * 180 / Pival, p4_Km[i].P());
                    if (p4_Km[i].P() > 0)
                        hist_all_Km_p_vs_phi->Fill(p4_Km[i].Phi() * 180 / Pival, p4_Km[i].P());
                    if (p4_Km[i].Theta() > 0 && p4_Km[i].Phi() != 0)
                        hist_all_Km_theta_vs_phi->Fill(p4_Km[i].Phi() * 180 / Pival, p4_Km[i].Theta() * 180 / Pival);

                    if (p4_phot[i].P() > 0)
                        hist_all_photon_p->Fill(p4_phot[i].P());
                    if (p4_phot[i].Theta() > 0)
                        hist_all_photon_theta->Fill(p4_phot[i].Theta() * 180 / Pival);
                    if (p4_phot[i].Phi() != 0)
                        hist_all_photon_phi->Fill(p4_phot[i].Phi() * 180 / Pival);
                    if (p4_phot[i].P() > 0)
                        hist_all_photon_p_vs_theta->Fill(p4_phot[i].Theta() * 180 / Pival, p4_phot[i].P());
                    if (p4_phot[i].P() > 0)
                        hist_all_photon_p_vs_phi->Fill(p4_phot[i].Phi() * 180 / Pival, p4_phot[i].P());
                    if (p4_phot[i].Theta() > 0 && p4_phot[i].Phi() != 0)
                        hist_all_photon_theta_vs_phi->Fill(p4_phot[i].Phi() * 180 / Pival, p4_phot[i].Theta() * 180 / Pival);
                }

                for (Int_t i = 0; i < BUFFER; i++)
                {
                    if (partLUND_pid[i] != 0)
                        hist_LUND_pid->Fill(partLUND_pid[i]);
                }

                /// ///////////////////////////////////////////////////
                /// count generated particles of each type

                for (Int_t i = 0; i < BUFFER; i++)
                {
                    if (partLUND_pid[i] == 11)
                    {
                        e_MCcount += 1;
                    }
                    if (partLUND_pid[i] == 2212)
                    {
                        p_MCcount += 1;
                    }
                    if (partLUND_pid[i] == 2112)
                    {
                        n_MCcount += 1;
                    }
                    if (partLUND_pid[i] == 211)
                    {
                        pip_MCcount += 1;
                    }
                    if (partLUND_pid[i] == -211)
                    {
                        pim_MCcount += 1;
                    }
                    if (partLUND_pid[i] == 321)
                    {
                        Kp_MCcount += 1;
                    }
                    if (partLUND_pid[i] == -321)
                    {
                        Km_MCcount += 1;
                    }
                }

                int e_RECcount = 0;
                int p_RECcount = 0;
                int n_RECcount = 0;
                int pip_RECcount = 0;
                int pim_RECcount = 0;
                int Kp_RECcount = 0;
                int Km_RECcount = 0;

                for (Int_t i = 0; i < BUFFER; i++)
                {
                    if (part_pid[i] == 11)
                    {
                        e_RECcount += 1;
                    }
                    if (part_pid[i] == 2212)
                    {
                        p_RECcount += 1;
                    }
                    if (part_pid[i] == 2112)
                    {
                        n_RECcount += 1;
                    }
                    if (part_pid[i] == 211)
                    {
                        pip_RECcount += 1;
                    }
                    if (part_pid[i] == -211)
                    {
                        pim_RECcount += 1;
                    }
                    if (part_pid[i] == 321)
                    {
                        Kp_RECcount += 1;
                    }
                    if (part_pid[i] == -321)
                    {
                        Km_RECcount += 1;
                    }
                }


                /// //////////////////////////////////////////////////////////////////////
                /// track quality:

                for (Int_t i = 0; i < BUFFER; i++)
                {

                    if (part_DC_Track_chi2[i] > 0)
                        hist_FD_chi2->Fill(part_DC_Track_chi2[i]);
                    if (part_DC_Track_chi2[i] > 0)
                        hist_FD_chi2_NDF->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);

                    if (part_DC_sector[i] == 1)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[0]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[0]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                    if (part_DC_sector[i] == 2)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[1]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[1]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                    if (part_DC_sector[i] == 3)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[2]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[2]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                    if (part_DC_sector[i] == 4)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[3]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[3]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                    if (part_DC_sector[i] == 5)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[4]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[4]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                    if (part_DC_sector[i] == 6)
                    {
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_sec[5]->Fill(part_DC_Track_chi2[i]);
                        if (part_DC_Track_chi2[i] > 0)
                            hist_FD_chi2_NDF_sec[5]->Fill(part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                    }
                }

                /// ///////////////////////////////////////////////////////
                /// Calorimeter energy depositions
                ///
                /// e_FTOF_sec[i] e_PCAL_sec[i] pip_FTOF_sec[i] pim_FTOF_sec[i]

                for (Int_t i = 0; i < BUFFER; i++)
                {

                    if (part_Cal_energy_total[i] > 0 && FD_eid_all_check[i])
                        hist_CAL_Edep_electron->Fill(part_Cal_energy_total[i]);
                    if (part_Cal_energy_total[i] > 0 && FD_protid_all_check[i])
                        hist_CAL_Edep_proton->Fill(part_Cal_energy_total[i]);
                    if (part_Cal_energy_total[i] > 0 && FD_pipid_all_check[i])
                        hist_CAL_Edep_pip->Fill(part_Cal_energy_total[i]);
                    if (part_Cal_energy_total[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]))
                        hist_CAL_Edep_pim->Fill(part_Cal_energy_total[i]);
                    if (part_Cal_energy_total[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i]))
                        hist_CAL_Edep_pip_diffsecele->Fill(part_Cal_energy_total[i]);
                    if (part_Cal_energy_total[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i]))
                        hist_CAL_Edep_pim_diffsecele->Fill(part_Cal_energy_total[i]);

                    if (part_Cal_PCAL_energy[i] > 0 && FD_eid_all_check[i])
                        hist_PCAL_Edep_electron->Fill(part_Cal_PCAL_energy[i]);
                    if (part_Cal_PCAL_energy[i] > 0 && FD_protid_all_check[i])
                        hist_PCAL_Edep_proton->Fill(part_Cal_PCAL_energy[i]);
                    if (part_Cal_PCAL_energy[i] > 0 && FD_pipid_all_check[i])
                        hist_PCAL_Edep_pip->Fill(part_Cal_PCAL_energy[i]);
                    if (part_Cal_PCAL_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]))
                        hist_PCAL_Edep_pim->Fill(part_Cal_PCAL_energy[i]);
                    if (part_Cal_PCAL_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i]))
                        hist_PCAL_Edep_pip_diffsecele->Fill(part_Cal_PCAL_energy[i]);
                    if (part_Cal_PCAL_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i]))
                        hist_PCAL_Edep_pim_diffsecele->Fill(part_Cal_PCAL_energy[i]);

                    if (part_Cal_ECin_energy[i] > 0 && FD_eid_all_check[i])
                        hist_ECin_Edep_electron->Fill(part_Cal_ECin_energy[i]);
                    if (part_Cal_ECin_energy[i] > 0 && FD_protid_all_check[i])
                        hist_ECin_Edep_proton->Fill(part_Cal_ECin_energy[i]);
                    if (part_Cal_ECin_energy[i] > 0 && FD_pipid_all_check[i])
                        hist_ECin_Edep_pip->Fill(part_Cal_ECin_energy[i]);
                    if (part_Cal_ECin_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]))
                        hist_ECin_Edep_pim->Fill(part_Cal_ECin_energy[i]);
                    if (part_Cal_ECin_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i]))
                        hist_ECin_Edep_pip_diffsecele->Fill(part_Cal_ECin_energy[i]);
                    if (part_Cal_ECin_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i]))
                        hist_ECin_Edep_pim_diffsecele->Fill(part_Cal_ECin_energy[i]);

                    if (part_Cal_ECout_energy[i] > 0 && FD_eid_all_check[i])
                        hist_ECout_Edep_electron->Fill(part_Cal_ECout_energy[i]);
                    if (part_Cal_ECout_energy[i] > 0 && FD_protid_all_check[i])
                        hist_ECout_Edep_proton->Fill(part_Cal_ECout_energy[i]);
                    if (part_Cal_ECout_energy[i] > 0 && FD_pipid_all_check[i])
                        hist_ECout_Edep_pip->Fill(part_Cal_ECout_energy[i]);
                    if (part_Cal_ECout_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]))
                        hist_ECout_Edep_pim->Fill(part_Cal_ECout_energy[i]);
                    if (part_Cal_ECout_energy[i] > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i]))
                        hist_ECout_Edep_pip_diffsecele->Fill(part_Cal_ECout_energy[i]);
                    if (part_Cal_ECout_energy[i] > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i]))
                        hist_ECout_Edep_pim_diffsecele->Fill(part_Cal_ECout_energy[i]);

                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && FD_eid_all_check[i])
                        hist_ECAL_Edep_electron->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && FD_protid_all_check[i])
                        hist_ECAL_Edep_proton->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && FD_pipid_all_check[i])
                        hist_ECAL_Edep_pip->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]))
                        hist_ECAL_Edep_pim->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && FD_pipid_all_check[i] && (e_FTOF_sec[i] != pip_FTOF_sec[i]))
                        hist_ECAL_Edep_pip_diffsecele->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                    if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && (FD_pimid_default_PID_check[i] && FD_pimid_beta_check[i]) && (e_FTOF_sec[i] != pim_FTOF_sec[i]))
                        hist_ECAL_Edep_pim_diffsecele->Fill(part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                }


                /// ////////////////////////////////////////////////////////

                // fill vertex plots

                for (Int_t i = 0; i < BUFFER; i++)
                {

                    if (part_charge[i] == +1 && part_vz[i] != 0)
                        hist_positive_vertex->Fill(part_vz[i]);
                    if (part_charge[i] == -1 && part_vz[i] != 0)
                        hist_negative_vertex->Fill(part_vz[i]);
                    if (part_charge[i] == +1 && part_vz[i] != 0 && part_theta[i] > 0)
                        hist_positive_vertex_vs_theta->Fill(part_theta[i] * 180 / Pival, part_vz[i]);
                    if (part_charge[i] == -1 && part_vz[i] != 0 && part_theta[i] > 0)
                        hist_negative_vertex_vs_theta->Fill(part_theta[i] * 180 / Pival, part_vz[i]);
                    if (part_charge[i] == +1 && part_vz[i] != 0 && part_phi[i] != 0)
                        hist_positive_vertex_vs_phi->Fill(part_phi[i] * 180 / Pival, part_vz[i]);
                    if (part_charge[i] == -1 && part_vz[i] != 0 && part_phi[i] != 0)
                        hist_negative_vertex_vs_phi->Fill(part_phi[i] * 180 / Pival, part_vz[i]);
                    if (part_charge[i] == +1 && part_vz[i] != 0 && part_p[i] > 0)
                        hist_positive_vertex_vs_p->Fill(part_p[i], part_vz[i]);
                    if (part_charge[i] == -1 && part_vz[i] != 0 && part_p[i] > 0)
                        hist_negative_vertex_vs_p->Fill(part_p[i], part_vz[i]);

                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == +1 && part_vz[i] != 0)
                        hist_positive_vertex_sec[i]->Fill(part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == -1 && part_vz[i] != 0)
                        hist_negative_vertex_sec[i]->Fill(part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == +1 && part_vz[i] != 0 && part_theta[i] > 0)
                        hist_positive_vertex_vs_theta_sec[i]->Fill(part_theta[i] * 180 / Pival, part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == -1 && part_vz[i] != 0 && part_theta[i] > 0)
                        hist_negative_vertex_vs_theta_sec[i]->Fill(part_theta[i] * 180 / Pival, part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == +1 && part_vz[i] != 0 && part_phi[i] != 0)
                        hist_positive_vertex_vs_phi_sec[i]->Fill(part_phi[i] * 180 / Pival, part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == -1 && part_vz[i] != 0 && part_phi[i] != 0)
                        hist_negative_vertex_vs_phi_sec[i]->Fill(part_phi[i] * 180 / Pival, part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == +1 && part_vz[i] != 0 && part_p[i] > 0)
                        hist_positive_vertex_vs_p_sec[i]->Fill(part_p[i], part_vz[i]);
                    if (part_Cal_PCAL_sector[i] == i + 1 && part_charge[i] == -1 && part_vz[i] != 0 && part_p[i] > 0)
                        hist_negative_vertex_vs_p_sec[i]->Fill(part_p[i], part_vz[i]);
                }

                for (Int_t i = 0; i < BUFFER; i++)
                {

                    if (e_vz[i] != 0)
                        hist_electron_vertex->Fill(e_vz[i]);
                    if (p_vz[i] != 0)
                        hist_proton_vertex->Fill(p_vz[i]);
                    if (e_vz[i] != 0 && p4_ele[i].Theta() > 0)
                        hist_electron_vertex_vs_theta->Fill(p4_ele[i].Theta() * 180 / Pival, e_vz[i]);
                    if (p_vz[i] != 0 && p4_prot[i].Theta() > 0)
                        hist_proton_vertex_vs_theta->Fill(p4_prot[i].Theta() * 180 / Pival, e_vz[i]);
                    if (e_vz[i] != 0 && p4_ele[i].Phi() != 0)
                        hist_electron_vertex_vs_phi->Fill(p4_ele[i].Phi() * 180 / Pival, e_vz[i]);
                    if (p_vz[i] != 0 && p4_prot[i].Phi() != 0)
                        hist_proton_vertex_vs_phi->Fill(p4_prot[i].Phi() * 180 / Pival, e_vz[i]);
                    if (e_vz[i] != 0 && p4_ele[i].P() > 0)
                        hist_electron_vertex_vs_p->Fill(p4_ele[i].P(), e_vz[i]);
                    if (p_vz[i] != 0 && p4_prot[i].P() > 0)
                        hist_proton_vertex_vs_p->Fill(p4_prot[i].P(), e_vz[i]);

                    if (e_PCAL_sec[i] == i + 1 && e_vz[i] != 0)
                        hist_electron_vertex_sec[i]->Fill(e_vz[i]);
                    if (p_PCAL_sec[i] == i + 1 && p_vz[i] != 0)
                        hist_proton_vertex_sec[i]->Fill(p_vz[i]);
                    if (e_PCAL_sec[i] == i + 1 && e_vz[i] != 0 && p4_ele[i].Theta() > 0)
                        hist_electron_vertex_vs_theta_sec[i]->Fill(p4_ele[i].Theta() * 180 / Pival, e_vz[i]);
                    if (p_PCAL_sec[i] == i + 1 && p_vz[i] != 0 != 0 && p4_prot[i].Theta() > 0)
                        hist_proton_vertex_vs_theta_sec[i]->Fill(p4_prot[i].Theta() * 180 / Pival, e_vz[i]);
                    if (e_PCAL_sec[i] == i + 1 && e_vz[i] != 0 && p4_ele[i].Phi() != 0)
                        hist_electron_vertex_vs_phi_sec[i]->Fill(p4_ele[i].Phi() * 180 / Pival, e_vz[i]);
                    if (p_PCAL_sec[i] == i + 1 && p_vz[i] != 0 && p4_prot[i].Phi() != 0)
                        hist_proton_vertex_vs_phi_sec[i]->Fill(p4_prot[i].Phi() * 180 / Pival, e_vz[i]);
                    if (e_PCAL_sec[i] == i + 1 && e_vz[i] != 0 && p4_ele[i].P() > 0)
                        hist_electron_vertex_vs_p_sec[i]->Fill(p4_ele[i].P(), e_vz[i]);
                    if (p_PCAL_sec[i] == i + 1 && p_vz[i] != 0 && p4_prot[i].P() > 0)
                        hist_proton_vertex_vs_p_sec[i]->Fill(p4_prot[i].P(), e_vz[i]);
                }


                // Phi angle versus sector

                for (Int_t i = 0; i < BUFFER; i++)
                {
                    if (part_charge[i] == +1 && part_phi[i] != 0)
                        hist_FTOF_phi_vs_sector_positives->Fill(part_FTOF_sector_layer2[i], part_phi[i] * 180 / Pival);
                    if (part_charge[i] == -1 && part_phi[i] != 0)
                        hist_FTOF_phi_vs_sector_negatives->Fill(part_FTOF_sector_layer2[i], part_phi[i] * 180 / Pival);
                }
            }

            /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///  fill the PID historgrams:
            /// ///////////////////////////////////////////////////////////////

            for (Int_t i = 0; i < BUFFER; i++){
            
            
              if (FD_eid_default_PID_check[i]){
                if(part_Cal_PCAL_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_Cal_PCAL_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_Cal_PCAL_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_Cal_PCAL_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_Cal_PCAL_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_Cal_PCAL_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_ele->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_ele->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_ele->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }
              
              if (FD_protid_default_PID_check[i]){
                if(part_DC_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_prot->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_prot->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_prot->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }
              
              if (FD_pipid_default_PID_check[i]){
                if(part_DC_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_pip->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_pip->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_pip->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }
              
              if (FD_pimid_default_PID_check[i]){
                if(part_DC_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_pim->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_pim->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_pim->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }
              
              if (FD_Kpid_default_PID_check[i]){
                if(part_DC_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_Kp->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_Kp->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_Kp->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }
              
              if (FD_Kmid_default_PID_check[i]){
                if(part_DC_sector[i] == 1){
                  hist_DC_edge_chi2_region1_sec1_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec1_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec1_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 2){
                  hist_DC_edge_chi2_region1_sec2_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec2_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec2_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 3){
                  hist_DC_edge_chi2_region1_sec3_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec3_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec3_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 4){
                  hist_DC_edge_chi2_region1_sec4_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec4_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec4_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 5){
                  hist_DC_edge_chi2_region1_sec5_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec5_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec5_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
                if(part_DC_sector[i] == 6){
                  hist_DC_edge_chi2_region1_sec6_Km->Fill(part_DC_edge1[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region2_sec6_Km->Fill(part_DC_edge2[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                  hist_DC_edge_chi2_region3_sec6_Km->Fill(part_DC_edge3[i], part_DC_Track_chi2[i] / part_DC_Track_NDF[i]);
                }
              }




                if (fill_electron_pid_histograms)
                {

                    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // all particles with negative charge

                    if (FD_eid_charge_check[i])
                    { // only particles which passed the default pid cut

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[0]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[0]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[0]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[0]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[0]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[0]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[0]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[0]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[0]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[0]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[0]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[0]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[0]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[0]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[0]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[0]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[0]->Fill(part_vz[i]);
                    }

                    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // default PID  (neg. charge is fulfilled for all)


                    if (FD_pipid_default_PID_check[i])
                    {
                      if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec1->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                      if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec2->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                      if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0 && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec3->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                      if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec4->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                      if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec5->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                      if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                        hist_PCAL_vs_ECin_sampling_fraction_pion_sec6->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                      }
                    }




                    if (FD_eid_default_PID_check[i])
                    {

                        double theta_PCAL = 180 / Pival * acos(part_Cal_PCAL_z[i] / sqrt(pow(part_Cal_PCAL_x[i], 2) + pow(part_Cal_PCAL_y[i], 2) + pow(part_Cal_PCAL_z[i], 2)));
                        double phi_PCAL_raw = 180 / Pival * atan2(part_Cal_PCAL_y[i] / sqrt(pow(part_Cal_PCAL_x[i], 2) + pow(part_Cal_PCAL_y[i], 2) + pow(part_Cal_PCAL_z[i], 2)), part_Cal_PCAL_x[i] / sqrt(pow(part_Cal_PCAL_x[i], 2) + pow(part_Cal_PCAL_y[i], 2) + pow(part_Cal_PCAL_z[i], 2)));
                        double phi_PCAL = 0;
                        if (part_Cal_PCAL_sector[i] == 1)
                            phi_PCAL = phi_PCAL_raw;
                        if (part_Cal_PCAL_sector[i] == 2)
                            phi_PCAL = phi_PCAL_raw - 60;
                        if (part_Cal_PCAL_sector[i] == 3)
                            phi_PCAL = phi_PCAL_raw - 120;
                        if (part_Cal_PCAL_sector[i] == 4 && phi_PCAL_raw > 0)
                            phi_PCAL = phi_PCAL_raw - 180;
                        if (part_Cal_PCAL_sector[i] == 4 && phi_PCAL_raw < 0)
                            phi_PCAL = phi_PCAL_raw + 180;
                        if (part_Cal_PCAL_sector[i] == 5)
                            phi_PCAL = phi_PCAL_raw + 120;
                        if (part_Cal_PCAL_sector[i] == 6)
                            phi_PCAL = phi_PCAL_raw + 60;

                        for (Int_t j = 0; j < 30; j++)
                        {
                            if (theta_PCAL > j + 5 && theta_PCAL <= j + 6)
                            {
                                if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec1[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec2[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec3[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec4[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec5[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi_sec6[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                                if (part_p[i] > 0)
                                    hist_electron_sampfrac_vs_phi[j]->Fill(phi_PCAL, part_Cal_energy_total[i] / part_p[i]);
                            }
                        }


                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec1->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec1->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec1_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec1_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec2->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec2->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec2_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec2_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0 && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec3->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec3->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec3_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec3_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec4->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec4->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec4_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec4_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec5->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec5->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec5_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec5_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0  && part_Cal_ECin_energy[i] > 0){
                          hist_PCAL_vs_ECin_sampling_fraction_sec6->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(FD_eid_EC_sampling_fraction_check[i]) hist_PCAL_vs_ECin_sampling_fraction_cut_sec6->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                          if(part_p[i] > 2 && part_p[i] < 3){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p23->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 3 && part_p[i] < 4){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p34->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 4 && part_p[i] < 5){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p45->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 5 && part_p[i] < 6){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p56->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 6 && part_p[i] < 7){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p67->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 7 && part_p[i] < 8){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p78->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 8 && part_p[i] < 9){ hist_PCAL_vs_ECin_sampling_fraction_sec6_p89->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                          if(part_p[i] > 9){                  hist_PCAL_vs_ECin_sampling_fraction_sec6_p9->Fill(part_Cal_ECin_energy[i] / part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);}
                        }



                        double theta_DCr1 = 180 / Pival * acos(part_DC_c1z[i] / sqrt(pow(part_DC_c1x[i], 2) + pow(part_DC_c1y[i], 2) + pow(part_DC_c1z[i], 2)));
                        double phi_DCr1_raw = 180 / Pival * atan2(part_DC_c1y[i] / sqrt(pow(part_DC_c1x[i], 2) + pow(part_DC_c1y[i], 2) + pow(part_DC_c1z[i], 2)), part_DC_c1x[i] / sqrt(pow(part_DC_c1x[i], 2) + pow(part_DC_c1y[i], 2) + pow(part_DC_c1z[i], 2)));
                        double phi_DCr1 = 0;
                        if (part_DC_sector[i] == 1)
                            phi_DCr1 = phi_DCr1_raw;
                        if (part_DC_sector[i] == 2)
                            phi_DCr1 = phi_DCr1_raw - 60;
                        if (part_DC_sector[i] == 3)
                            phi_DCr1 = phi_DCr1_raw - 120;
                        if (part_DC_sector[i] == 4 && phi_DCr1_raw > 0)
                            phi_DCr1 = phi_DCr1_raw - 180;
                        if (part_DC_sector[i] == 4 && phi_DCr1_raw < 0)
                            phi_DCr1 = phi_DCr1_raw + 180;
                        if (part_DC_sector[i] == 5)
                            phi_DCr1 = phi_DCr1_raw + 120;
                        if (part_DC_sector[i] == 6)
                            phi_DCr1 = phi_DCr1_raw + 60;

                        for (Int_t j = 0; j < 30; j++)
                        {
                            if (theta_PCAL > j + 5 && theta_PCAL <= j + 6)
                            {
                                if (part_DC_sector[i] == 1 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec1[j]->Fill(phi_DCr1);
                                if (part_DC_sector[i] == 2 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec2[j]->Fill(phi_DCr1);
                                if (part_DC_sector[i] == 3 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec3[j]->Fill(phi_DCr1);
                                if (part_DC_sector[i] == 4 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec4[j]->Fill(phi_DCr1);
                                if (part_DC_sector[i] == 5 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec5[j]->Fill(phi_DCr1);
                                if (part_DC_sector[i] == 6 && part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi_sec6[j]->Fill(phi_DCr1);
                                if (part_p[i] > 0)
                                    hist_DC_region1_electron_counts_vs_phi[j]->Fill(phi_DCr1);
                            }
                        }

                        double theta_DCr2 = 180 / Pival * acos(part_DC_c2z[i] / sqrt(pow(part_DC_c2x[i], 2) + pow(part_DC_c2y[i], 2) + pow(part_DC_c2z[i], 2)));
                        double phi_DCr2_raw = 180 / Pival * atan2(part_DC_c2y[i] / sqrt(pow(part_DC_c2x[i], 2) + pow(part_DC_c2y[i], 2) + pow(part_DC_c2z[i], 2)), part_DC_c2x[i] / sqrt(pow(part_DC_c2x[i], 2) + pow(part_DC_c2y[i], 2) + pow(part_DC_c2z[i], 2)));
                        double phi_DCr2 = 0;
                        if (part_DC_sector[i] == 1)
                            phi_DCr2 = phi_DCr2_raw;
                        if (part_DC_sector[i] == 2)
                            phi_DCr2 = phi_DCr2_raw - 60;
                        if (part_DC_sector[i] == 3)
                            phi_DCr2 = phi_DCr2_raw - 120;
                        if (part_DC_sector[i] == 4 && phi_DCr2_raw > 0)
                            phi_DCr2 = phi_DCr2_raw - 180;
                        if (part_DC_sector[i] == 4 && phi_DCr2_raw < 0)
                            phi_DCr2 = phi_DCr2_raw + 180;
                        if (part_DC_sector[i] == 5)
                            phi_DCr2 = phi_DCr2_raw + 120;
                        if (part_DC_sector[i] == 6)
                            phi_DCr2 = phi_DCr2_raw + 60;

                        for (Int_t j = 0; j < 30; j++)
                        {
                            if (theta_PCAL > j + 5 && theta_PCAL <= j + 6)
                            {
                                if (part_DC_sector[i] == 1 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec1[j]->Fill(phi_DCr2);
                                if (part_DC_sector[i] == 2 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec2[j]->Fill(phi_DCr2);
                                if (part_DC_sector[i] == 3 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec3[j]->Fill(phi_DCr2);
                                if (part_DC_sector[i] == 4 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec4[j]->Fill(phi_DCr2);
                                if (part_DC_sector[i] == 5 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec5[j]->Fill(phi_DCr2);
                                if (part_DC_sector[i] == 6 && part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi_sec6[j]->Fill(phi_DCr2);
                                if (part_p[i] > 0)
                                    hist_DC_region2_electron_counts_vs_phi[j]->Fill(phi_DCr2);
                            }
                        }

                        double theta_DCr3 = 180 / Pival * acos(part_DC_c3z[i] / sqrt(pow(part_DC_c3x[i], 2) + pow(part_DC_c3y[i], 2) + pow(part_DC_c3z[i], 2)));
                        double phi_DCr3_raw = 180 / Pival * atan2(part_DC_c3y[i] / sqrt(pow(part_DC_c3x[i], 2) + pow(part_DC_c3y[i], 2) + pow(part_DC_c3z[i], 2)), part_DC_c3x[i] / sqrt(pow(part_DC_c3x[i], 2) + pow(part_DC_c3y[i], 2) + pow(part_DC_c3z[i], 2)));
                        double phi_DCr3 = 0;
                        if (part_DC_sector[i] == 1)
                            phi_DCr3 = phi_DCr3_raw;
                        if (part_DC_sector[i] == 2)
                            phi_DCr3 = phi_DCr3_raw - 60;
                        if (part_DC_sector[i] == 3)
                            phi_DCr3 = phi_DCr3_raw - 120;
                        if (part_DC_sector[i] == 4 && phi_DCr3_raw > 0)
                            phi_DCr3 = phi_DCr3_raw - 180;
                        if (part_DC_sector[i] == 4 && phi_DCr3_raw < 0)
                            phi_DCr3 = phi_DCr3_raw + 180;
                        if (part_DC_sector[i] == 5)
                            phi_DCr3 = phi_DCr3_raw + 120;
                        if (part_DC_sector[i] == 6)
                            phi_DCr3 = phi_DCr3_raw + 60;

                        for (Int_t j = 0; j < 30; j++)
                        {
                            if (theta_PCAL > j + 5 && theta_PCAL <= j + 6)
                            {
                                if (part_DC_sector[i] == 1 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec1[j]->Fill(phi_DCr3);
                                if (part_DC_sector[i] == 2 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec2[j]->Fill(phi_DCr3);
                                if (part_DC_sector[i] == 3 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec3[j]->Fill(phi_DCr3);
                                if (part_DC_sector[i] == 4 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec4[j]->Fill(phi_DCr3);
                                if (part_DC_sector[i] == 5 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec5[j]->Fill(phi_DCr3);
                                if (part_DC_sector[i] == 6 && part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi_sec6[j]->Fill(phi_DCr3);
                                if (part_p[i] > 0)
                                    hist_DC_region3_electron_counts_vs_phi[j]->Fill(phi_DCr3);
                            }
                        }

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[1]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[1]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[1]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[1]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[1]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[1]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[1]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[1]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[1]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[1]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[1]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[1]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[1]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[1]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[1]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[1]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[1]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[1]->Fill(part_vz[i]);
                    }

                    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // EC cuts

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[2]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[2]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[2]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[2]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[2]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[2]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[2]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[2]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[2]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[2]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[2]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[2]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[2]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[2]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[2]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[2]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[2]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[2]->Fill(part_vz[i]);
                    }

                    // if(FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_DC_z_vertex_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_EC_sampling_fraction_check[i]){

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_sampling_fraction_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[3]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[3]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[3]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[3]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[3]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[3]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[3]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[3]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[3]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[3]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[3]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[3]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[3]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[3]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[3]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[3]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[3]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[3]->Fill(part_vz[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_hit_position_fiducial_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[4]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[4]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[4]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[4]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[4]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[4]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[4]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[4]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[4]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[4]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[4]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[4]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[4]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[4]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[4]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[4]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[4]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[4]->Fill(part_vz[i]);
                    }

                    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // DC cuts

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && DC_fiducial_cut_edge(i, 1) == true)
                    {
                        hist_DC_hit_position_region1_new_cut->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                    }
                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && DC_fiducial_cut_edge(i, 2) == true)
                    {
                        hist_DC_hit_position_region2_new_cut->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && DC_fiducial_cut_edge(i, 3) == true)
                    {
                        hist_DC_hit_position_region3_new_cut->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[5]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[5]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[5]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[5]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[5]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[5]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[5]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[5]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[5]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[5]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[5]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[5]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[5]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[5]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[5]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[5]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[5]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[5]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[5]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[5]->Fill(part_vz[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i])
                    {
                        hist_DC_hit_position_region2_cut5a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[6]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[6]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[6]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[6]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[6]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[6]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[6]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[6]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[6]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[6]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[6]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[6]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[6]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[6]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[6]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[6]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[6]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[6]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[6]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[6]->Fill(part_vz[i]);
                    }

                    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // vertex cut

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_DC_z_vertex_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[7]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[7]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[7]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[7]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[7]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[7]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[7]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[7]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[7]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[7]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[7]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[7]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[7]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[7]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[7]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[7]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[7]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[7]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[7]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[7]->Fill(part_vz[i]);
                    }

                    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // all other

                    // CC

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {
                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[8]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[8]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[8]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);
                    }

                    // EC

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {

                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[8]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[8]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {
                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[8]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[8]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[8]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[8]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[8]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[8]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);
                    }

                    // DC

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {

                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[8]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {

                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                    {

                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[8]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                    }

                    if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i])
                    {
                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[8]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[8]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[8]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[8]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[8]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[8]->Fill(part_vz[i]);
                    }

                    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // all cuts

                    if (FD_eid_all_check[i])
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[9]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[9]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[9]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[9]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[9]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[9]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[9]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[9]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[9]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[9]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[9]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[9]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[9]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[9]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[9]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[9]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[9]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[9]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[9]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[9]->Fill(part_vz[i]);

                    }

                    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // inverse cut  (neg. charge but no PID 11)

                    if (FD_eid_default_PID_check[i] == false && FD_eid_charge_check[i] == true)
                    {

                        if (part_p[i] > 0)
                            hist_HTCC_theta_vs_phi[10]->Fill(part_CC_HTCC_phi[i] * 180 / Pival, part_CC_HTCC_theta[i] * 180 / Pival);
                        if (part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_nphe[10]->Fill(part_CC_HTCC_nphe[i]);
                        if (part_CC_HTCC_nphe[i] > 0 && part_p[i] > 0 && part_Cal_energy_total[i] / part_p[i] > 0)
                            hist_HTCC_nphe_vs_sampling_fraction[10]->Fill(part_CC_HTCC_nphe[i], part_Cal_energy_total[i] / part_p[i]);

                        // EC
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL[10]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_ECin_energy[i] > 0 && part_Cal_ECout_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner[10]->Fill(part_Cal_ECin_energy[i], part_Cal_ECout_energy[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec1[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec2[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec3[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec4[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec5[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_total_sampling_fraction_sec6[10]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec1[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec2[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec3[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec4[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec5[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_PCAL_sampling_fraction_sec6[10]->Fill(part_p[i], part_Cal_PCAL_energy[i] / part_p[i]);

                        if (part_Cal_ECin_sector[i] == 1 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec1[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 2 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec2[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 3 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec3[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 4 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec4[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 5 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec5[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);
                        if (part_Cal_ECin_sector[i] == 6 && part_p[i] > 0)
                            hist_EC_ECAL_sampling_fraction_sec6[10]->Fill(part_p[i], (part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) / part_p[i]);

                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position[10]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                        if (part_Cal_ECin_x[i] != 0 && part_Cal_ECin_y[i] != 0)
                            hist_EC_inner_hit_position[10]->Fill(part_Cal_ECin_x[i], part_Cal_ECin_y[i]);
                        if (part_Cal_ECout_x[i] != 0 && part_Cal_ECout_y[i] != 0)
                            hist_EC_outer_hit_position[10]->Fill(part_Cal_ECout_x[i], part_Cal_ECout_y[i]);

                        // DC
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1[10]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2[10]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3[10]->Fill(part_DC_c3x[i], part_DC_c3y[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec1[10]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec2[10]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec3[10]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec4[10]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec5[10]->Fill(part_vz[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_vz[i] != 0)
                            hist_DC_z_vertex_sec6[10]->Fill(part_vz[i]);
                    }

                } // end of fill electron histogram selection

                /// //////////////////////////////////////////////////////////////////////////////////////////////////////
                /// TOF charged hadron plots
                /// /////////////////////////////////////////////////////////////////

                double beta_charge = Beta_charged(i, run);
                double beta_neutr = Beta_neutral(i, run);

                if (fill_hadron_pid_histograms)
                {

                    // proton

                    if (FD_protid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[0]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[0]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[0]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[0]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[0]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[0]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[0]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[0]->Fill(Getdvz(i));
                    }
                    if (FD_protid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[1]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[1]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[1]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0.4 && beta_charge > 0)
                            hist_beta_vs_p[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[1]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[1]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[1]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[1]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[1]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[2]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[2]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[2]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[2]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[2]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[2]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[2]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[2]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i])
                    {
                        hist_DC_hit_position_region2_hadron_cut_02a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[3]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[3]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[3]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[3]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[3]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[3]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[3]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[3]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[4]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[4]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[4]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[4]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[4]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[4]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[4]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[4]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[5]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[5]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[5]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[5]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[5]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[5]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[5]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[5]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[5]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[6]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[6]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[6]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[6]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[6]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[6]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[6]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[6]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[6]->Fill(Getdvz(i));
                    }
                    if (FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_maximum_probability_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[7]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[7]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[7]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[7]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[7]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[7]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[7]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[7]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[7]->Fill(Getdvz(i));
                    }
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[8]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[8]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[8]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[8]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[8]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[8]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[8]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[8]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[8]->Fill(Getdvz(i));
                    }
                    if (FD_protid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[9]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[9]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[9]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[9]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[9]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[9]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[9]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[9]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_p * m_p) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[9]->Fill(Getdvz(i));
                    }

                    // Neutron

                    if (FD_neutrid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[10]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[10]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[10]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[10]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[10]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[10]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[10]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[10]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[10]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[11]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[11]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[11]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[11]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[11]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[11]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[11]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[11]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[11]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[12]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[12]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[12]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[12]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[12]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[12]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[12]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[12]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[12]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[13]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[13]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[13]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[13]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[13]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[13]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[13]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[13]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[13]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[14]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[14]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[14]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[14]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[14]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[14]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[14]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[14]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[14]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[15]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[15]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[15]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[15]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[15]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[15]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[15]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[15]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[15]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[16]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[16]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[16]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[16]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[16]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[16]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[16]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[16]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[16]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[18]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[18]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[18]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[18]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[18]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[18]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[18]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[18]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[18]->Fill(Getdvz(i));
                    }
                    if (FD_neutrid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[19]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[19]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[19]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[19]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec1[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec2[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec3[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec4[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec5[19]->Fill(part_p[i], beta_neutr);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_sec6[19]->Fill(part_p[i], beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta_vs_p[19]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_delta_beta[19]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_n * m_n) - beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[19]->Fill(Getdvz(i));
                    }

                    // pip

                    if (FD_pipid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[20]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[20]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[20]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[20]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[20]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[20]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[20]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[20]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[20]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[21]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[21]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[21]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[21]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[21]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[21]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[21]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[21]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[21]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[22]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[22]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[22]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[22]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[22]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[22]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[22]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[22]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[22]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i])
                    {
                        hist_DC_hit_position_region2_hadron_cut_22a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[23]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[23]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[23]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[23]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[23]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[23]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[23]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[23]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[23]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[24]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[24]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[24]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[24]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[24]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[24]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[24]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[24]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[24]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[25]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[25]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[25]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[25]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[25]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[25]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[25]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[25]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[25]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[26]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[26]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[26]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[26]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[26]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[26]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[26]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[26]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[26]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_maximum_probability_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[27]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[27]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[27]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[27]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[27]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[27]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[27]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[27]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[27]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[28]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[28]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[28]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[28]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[28]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[28]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[28]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[28]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[28]->Fill(Getdvz(i));
                    }
                    if (FD_pipid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[29]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[29]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[29]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[29]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[29]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[29]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[29]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[29]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pip * m_pip) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[29]->Fill(Getdvz(i));
                    }

                    // pim

                    if (FD_pimid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[30]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[30]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[30]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[30]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[30]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[30]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[30]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[30]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[30]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_charge_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[31]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[31]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[31]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[31]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0.3 && beta_charge > 0)
                            hist_beta_vs_p[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[31]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[31]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[31]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[31]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[31]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[32]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[32]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[32]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[32]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[32]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[32]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[32]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[32]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[32]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i])
                    {
                        hist_DC_hit_position_region2_hadron_cut_32a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[33]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[33]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[33]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[33]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[33]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[33]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[33]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[33]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[33]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[34]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[34]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[34]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[34]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[34]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[34]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[34]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[34]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[34]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[35]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[35]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[35]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[35]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[35]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[35]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[35]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[35]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[35]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[36]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[36]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[36]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[36]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[36]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[36]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[36]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[36]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[36]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pipid_maximum_probability_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[37]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[37]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[37]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[37]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[37]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[37]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[37]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[37]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[37]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_ele_reject_check[i] && FD_pimid_EC_outer_vs_EC_inner_check[i] && FD_pimid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[38]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[38]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[38]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[38]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[38]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[38]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[38]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[38]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[38]->Fill(Getdvz(i));
                    }
                    if (FD_pimid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[39]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[39]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[39]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[39]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[39]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[39]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[39]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[39]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_pim * m_pim) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[39]->Fill(Getdvz(i));
                    }

                    // Kp

                    if (FD_Kpid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[40]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[40]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[40]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[40]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[40]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[40]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[40]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[40]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[40]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_charge_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[41]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[41]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[41]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[41]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[41]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[41]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[41]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[41]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[41]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[42]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[42]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[42]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[42]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[42]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[42]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[42]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[42]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[42]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i])
                    {
                        hist_DC_hit_position_region2_hadron_cut_42a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[43]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[43]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[43]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[43]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[43]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[43]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[43]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[43]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[43]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[44]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[44]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[44]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[44]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[44]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[44]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[44]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[44]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[44]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[45]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[45]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[45]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[45]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[45]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[45]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[45]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[45]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[45]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[46]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[46]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[46]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[46]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[46]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[46]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[46]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[46]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[46]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_maximum_probability_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[47]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[47]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[47]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[47]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[47]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[47]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[47]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[47]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[47]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[48]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[48]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[48]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[48]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[48]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[48]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[48]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[48]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[48]->Fill(Getdvz(i));
                    }
                    if (FD_Kpid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[49]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[49]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[49]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[49]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[49]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[49]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[49]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[49]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Kp * m_Kp) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[49]->Fill(Getdvz(i));
                    }

                    // Km

                    if (FD_Kmid_default_PID_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[50]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[50]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[50]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[50]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[50]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[50]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[50]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[50]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[50]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[51]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[51]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[51]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[51]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[51]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[51]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[51]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[51]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[51]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[52]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[52]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[52]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[52]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[52]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[52]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[52]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[52]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[52]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i])
                    {
                        hist_DC_hit_position_region2_hadron_cut_52a->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[53]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[53]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[53]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[53]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[53]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[53]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[53]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[53]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[53]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[54]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[54]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[54]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[54]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[54]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[54]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[54]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[54]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[54]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_delta_beta_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[55]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[55]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[55]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[55]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[55]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[55]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[55]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[55]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[55]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_tofmass_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[56]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[56]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[56]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[56]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[56]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[56]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[56]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[56]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[56]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_maximum_probability_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[57]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[57]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[57]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[57]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[57]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[57]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[57]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[57]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[57]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_ele_reject_check[i] && FD_Kmid_EC_outer_vs_EC_inner_check[i] && FD_Kmid_delta_vz_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[58]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[58]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[58]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[58]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[58]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[58]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[58]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[58]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[58]->Fill(Getdvz(i));
                    }
                    if (FD_Kmid_all_check[i])
                    {
                        if (part_DC_c1x[i] != 0 && part_DC_c1y[i] != 0)
                            hist_DC_hit_position_region1_hadron[59]->Fill(part_DC_c1x[i], part_DC_c1y[i]);
                        if (part_DC_c2x[i] != 0 && part_DC_c2y[i] != 0)
                            hist_DC_hit_position_region2_hadron[59]->Fill(part_DC_c2x[i], part_DC_c2y[i]);
                        if (part_DC_c3x[i] != 0 && part_DC_c3y[i] != 0)
                            hist_DC_hit_position_region3_hadron[59]->Fill(part_DC_c3x[i], part_DC_c3y[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_outer_vs_EC_inner_hadron[59]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec1[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec2[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec3[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec4[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec5[59]->Fill(part_p[i], beta_charge);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && beta_charge > 0)
                            hist_beta_vs_p_sec6[59]->Fill(part_p[i], beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta_vs_p[59]->Fill(part_p[i], part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (part_p[i] > 0 && beta_charge > 0)
                            hist_delta_beta[59]->Fill(part_p[i] / sqrt(part_p[i] * part_p[i] + m_Km * m_Km) - beta_charge);
                        if (Getdvz(i) != 0)
                            hist_delta_vz[59]->Fill(Getdvz(i));
                    }


                    /// /////////////////////////////////////////////////////////////////////////////////////////
                    /// Central detector

                    double beta_charge_central = Beta_charged_central(i, run);

                    // proton

                    if (CD_protid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[0]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[0]->Fill(Getdvz(i));
                    }
                    if (CD_protid_charge_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[1]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[1]->Fill(Getdvz(i));
                    }
                    if (CD_protid_charge_check[i] && CD_protid_beta_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[2]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[2]->Fill(Getdvz(i));
                    }
                    if (CD_protid_charge_check[i] && CD_protid_maximum_probability_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[3]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[3]->Fill(Getdvz(i));
                    }
                    if (CD_protid_charge_check[i] && CD_protid_delta_vz_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[4]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[4]->Fill(Getdvz(i));
                    }
                    if (CD_protid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[5]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[5]->Fill(Getdvz(i));
                    }

                    // Neutron

                    if (CD_neutrid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_CD_beta_vs_p[10]->Fill(part_p[i], beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[10]->Fill(Getdvz(i));
                    }
                    if (CD_neutrid_charge_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_CD_beta_vs_p[11]->Fill(part_p[i], beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[11]->Fill(Getdvz(i));
                    }
                    if (CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_beta_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_CD_beta_vs_p[12]->Fill(part_p[i], beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[12]->Fill(Getdvz(i));
                    }
                    if (CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_delta_vz_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_CD_beta_vs_p[14]->Fill(part_p[i], beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[14]->Fill(Getdvz(i));
                    }
                    if (CD_neutrid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_CD_beta_vs_p[15]->Fill(part_p[i], beta_neutr);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[15]->Fill(Getdvz(i));
                    }

                    // pip

                    if (CD_pipid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[20]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[20]->Fill(Getdvz(i));
                    }
                    if (CD_pipid_charge_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[21]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[21]->Fill(Getdvz(i));
                    }
                    if (CD_pipid_charge_check[i] && CD_pipid_beta_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[22]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[22]->Fill(Getdvz(i));
                    }
                    if (CD_pipid_charge_check[i] && CD_pipid_maximum_probability_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[23]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[23]->Fill(Getdvz(i));
                    }
                    if (CD_pipid_charge_check[i] && CD_pipid_delta_vz_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[24]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[24]->Fill(Getdvz(i));
                    }
                    if (CD_pipid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[25]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[25]->Fill(Getdvz(i));
                    }

                    // pim

                    if (CD_pimid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[30]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[30]->Fill(Getdvz(i));
                    }
                    if (CD_pimid_charge_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[31]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[31]->Fill(Getdvz(i));
                    }
                    if (CD_pimid_charge_check[i] && CD_pimid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[32]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[32]->Fill(Getdvz(i));
                    }
                    if (CD_pimid_charge_check[i] && CD_pipid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[33]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[33]->Fill(Getdvz(i));
                    }
                    if (CD_pimid_charge_check[i] && CD_pimid_delta_vz_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[34]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[34]->Fill(Getdvz(i));
                    }
                    if (CD_pimid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[35]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[35]->Fill(Getdvz(i));
                    }

                    // Kp

                    if (CD_Kpid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[40]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[40]->Fill(Getdvz(i));
                    }
                    if (CD_Kpid_charge_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[41]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[41]->Fill(Getdvz(i));
                    }
                    if (CD_Kpid_charge_check[i] && CD_Kpid_beta_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[42]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[42]->Fill(Getdvz(i));
                    }
                    if (CD_Kpid_charge_check[i] && CD_Kpid_maximum_probability_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[43]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[43]->Fill(Getdvz(i));
                    }
                    if (CD_Kpid_charge_check[i] && CD_Kpid_delta_vz_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[44]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[44]->Fill(Getdvz(i));
                    }
                    if (CD_Kpid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[45]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[45]->Fill(Getdvz(i));
                    }

                    // Km

                    if (CD_Kmid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[50]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[50]->Fill(Getdvz(i));
                    }
                    if (CD_Kmid_charge_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[51]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[51]->Fill(Getdvz(i));
                    }
                    if (CD_Kmid_charge_check[i] && CD_Kmid_beta_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[52]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[52]->Fill(Getdvz(i));
                    }
                    if (CD_Kmid_charge_check[i] && CD_Kmid_maximum_probability_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[53]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[53]->Fill(Getdvz(i));
                    }
                    if (CD_Kmid_charge_check[i] && CD_Kmid_delta_vz_check[i] && (beta_charge_central < 0.9999 || beta_charge_central > 1.0001))
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[54]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[54]->Fill(Getdvz(i));
                    }
                    if (CD_Kmid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_charge_central > 0)
                            hist_CD_beta_vs_p[55]->Fill(part_p[i], beta_charge_central);
                        if (Getdvz(i) != 0)
                            hist_CD_delta_vz[55]->Fill(Getdvz(i));
                    }

                } // end of fill hadron pid histogram selection

                /// //////////////////////////////////////////////////////////////////////////////////////////////////////
                /// photon ID plots
                /// /////////////////////////////////////////////////////////////////

                if (fill_photon_pid_histograms)
                {

                    if (FD_photid_default_PID_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_phot[0]->Fill(part_p[i], beta_neutr);
                        if (beta_neutr > 0)
                            hist_beta_phot[0]->Fill(beta_neutr);
                        if (part_p[i] > 0)
                            hist_EC_sampling_fraction_phot[0]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL_phot[0]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position_phot[0]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                    }

                    if (FD_photid_charge_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_phot[1]->Fill(part_p[i], beta_neutr);
                        if (beta_neutr > 0)
                            hist_beta_phot[1]->Fill(beta_neutr);
                        if (part_p[i] > 0)
                            hist_EC_sampling_fraction_phot[1]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL_phot[1]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position_phot[1]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                    }

                    if (FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_beta_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_phot[2]->Fill(part_p[i], beta_neutr);
                        if (beta_neutr > 0)
                            hist_beta_phot[2]->Fill(beta_neutr);
                        if (part_p[i] > 0)
                            hist_EC_sampling_fraction_phot[2]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL_phot[2]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position_phot[2]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                    }

                    if (FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_EC_hit_position_fiducial_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_phot[3]->Fill(part_p[i], beta_neutr);
                        if (beta_neutr > 0)
                            hist_beta_phot[3]->Fill(beta_neutr);
                        if (part_p[i] > 0)
                            hist_EC_sampling_fraction_phot[3]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL_phot[3]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position_phot[3]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                    }

                    if (FD_photid_all_check[i])
                    {
                        if (part_p[i] > 0 && beta_neutr > 0)
                            hist_beta_vs_p_phot[4]->Fill(part_p[i], beta_neutr);
                        if (beta_neutr > 0)
                            hist_beta_phot[4]->Fill(beta_neutr);
                        if (part_p[i] > 0)
                            hist_EC_sampling_fraction_phot[4]->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                        if ((part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]) > 0 && part_Cal_PCAL_energy[i] > 0)
                            hist_EC_PCAL_vs_EC_ECAL_phot[4]->Fill(part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i]);
                        if (part_Cal_PCAL_x[i] != 0 && part_Cal_PCAL_y[i] != 0)
                            hist_EC_PCAL_hit_position_phot[4]->Fill(part_Cal_PCAL_x[i], part_Cal_PCAL_y[i]);
                    }

                } // end of fill photon ID histogram selection

                /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /// Fill the FT histograms:

                // electrons

                if (FT_eid_charge_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[0]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[0]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[0]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[0]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[0]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_PID_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[1]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[1]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[1]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[1]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[1]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_FTCAL_fiducial_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[2]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[2]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[2]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[2]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[2]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_FTTRK_fiducial_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[3]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[3]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[3]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[3]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[3]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_FTHODO_fiducial_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[4]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[4]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[4]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[4]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[4]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_energy_vs_radius_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[5]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[5]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[5]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[5]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[5]->Fill(Beta_charged_FT(i, run));
                }

                if (FT_eid_all_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[6]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[6]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[6]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[6]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[6]->Fill(Beta_charged_FT(i, run));
                }

                hist_FT_FTCAL_energy_vs_radius[7]->Fill(part_FT_energy[i], part_FT_radius[i]);
                hist_FT_FTCAL_hit_position[7]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                hist_FT_FTTRK_hit_position[7]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                hist_FT_FTHODO_hit_position[7]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                hist_FT_beta[7]->Fill(Beta_charged_FT(i, run));

                // photons

                if (FT_photid_charge_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[10]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[10]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[10]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[10]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[10]->Fill(Beta_neutral_FT(i, run));
                }

                if (FT_photid_PID_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[11]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[11]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[11]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[11]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[11]->Fill(Beta_neutral_FT(i, run));
                }

                if (FT_photid_FTCAL_fiducial_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[12]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[12]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[12]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[12]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[12]->Fill(Beta_neutral_FT(i, run));
                }

                if (FT_photid_beta_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[13]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[13]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[13]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[13]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[13]->Fill(Beta_neutral_FT(i, run));
                }

                if (FT_photid_all_check[i])
                {
                    hist_FT_FTCAL_energy_vs_radius[14]->Fill(part_FT_energy[i], part_FT_radius[i]);
                    hist_FT_FTCAL_hit_position[14]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                    hist_FT_FTTRK_hit_position[14]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                    hist_FT_FTHODO_hit_position[14]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                    hist_FT_beta[14]->Fill(Beta_neutral_FT(i, run));
                }

                hist_FT_FTCAL_energy_vs_radius[15]->Fill(part_FT_energy[i], part_FT_radius[i]);
                hist_FT_FTCAL_hit_position[15]->Fill(part_FTCAL_x[i], part_FTCAL_y[i]);
                hist_FT_FTTRK_hit_position[15]->Fill(part_FTTRK_x[i], part_FTTRK_y[i]);
                hist_FT_FTHODO_hit_position[15]->Fill(part_FTHODO_x[i], part_FTHODO_y[i]);
                hist_FT_beta[15]->Fill(Beta_neutral_FT(i, run));

                /// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                /// define calorimeter fiducial cut borders

                if (FD_photid_default_PID_check[i])
                {

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_photon_u_coord_sec1->Fill(part_Cal_PCAL_lu[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_photon_u_coord_sec2->Fill(part_Cal_PCAL_lu[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_photon_u_coord_sec3->Fill(part_Cal_PCAL_lu[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_photon_u_coord_sec4->Fill(part_Cal_PCAL_lu[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_photon_u_coord_sec5->Fill(part_Cal_PCAL_lu[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_photon_u_coord_sec6->Fill(part_Cal_PCAL_lu[i]);
                    if (part_p[i] > 0)
                        hist_photon_u_coord->Fill(part_Cal_PCAL_lu[i]);

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_photon_v_coord_sec1->Fill(part_Cal_PCAL_lv[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_photon_v_coord_sec2->Fill(part_Cal_PCAL_lv[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_photon_v_coord_sec3->Fill(part_Cal_PCAL_lv[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_photon_v_coord_sec4->Fill(part_Cal_PCAL_lv[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_photon_v_coord_sec5->Fill(part_Cal_PCAL_lv[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_photon_v_coord_sec6->Fill(part_Cal_PCAL_lv[i]);
                    if (part_p[i] > 0)
                        hist_photon_v_coord->Fill(part_Cal_PCAL_lv[i]);

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_photon_w_coord_sec1->Fill(part_Cal_PCAL_lw[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_photon_w_coord_sec2->Fill(part_Cal_PCAL_lw[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_photon_w_coord_sec3->Fill(part_Cal_PCAL_lw[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_photon_w_coord_sec4->Fill(part_Cal_PCAL_lw[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_photon_w_coord_sec5->Fill(part_Cal_PCAL_lw[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_photon_w_coord_sec6->Fill(part_Cal_PCAL_lw[i]);
                    if (part_p[i] > 0)
                        hist_photon_w_coord->Fill(part_Cal_PCAL_lw[i]);
                }

                if (FD_eid_default_PID_check[i])
                {

                    if (part_p[i] > 3 && part_p[i] < 6)
                    {

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec1->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec2->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec3->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec4->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec5->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord_sec6->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_p[i] > 0)
                            hist_electron_sampfrac_vs_u_coord->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec1->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec2->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec3->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec4->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec5->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord_sec6->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_p[i] > 0)
                            hist_electron_sampfrac_vs_v_coord->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i]);

                        if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec1->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec2->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec3->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec4->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec5->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord_sec6->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                        if (part_p[i] > 0)
                            hist_electron_sampfrac_vs_w_coord->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i]);
                    }

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec1->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec2->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec3->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec4->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec5->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_u_coord_vs_p_sec6->Fill(part_Cal_PCAL_lu[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec1->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec2->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec3->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec4->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec5->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_v_coord_vs_p_sec6->Fill(part_Cal_PCAL_lv[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec1->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec2->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec3->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec4->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec5->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_electron_sampfrac_vs_w_coord_vs_p_sec6->Fill(part_Cal_PCAL_lw[i], part_Cal_energy_total[i] / part_p[i], part_p[i]);
                }

                if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                    hist_HTCC_Nphe_vs_momentum_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);

                if (FD_protid_all_check[i])
                {
                    if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_prot_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                }

                if (FD_pipid_all_check[i])
                {
                    if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec1->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec2->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec3->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec4->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec5->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                    if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_momentum_pip_sec6->Fill(part_p[i], part_CC_HTCC_nphe[i]);
                }

                if (part_p[i] > 4)
                {

                    if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                        hist_HTCC_Nphe_vs_beta_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);

                    if (FD_protid_all_check[i])
                    {
                        if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_prot_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    }

                    if (FD_pipid_all_check[i])
                    {
                        if (part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 1 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec1->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 2 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec2->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 3 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec3->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 4 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec4->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 5 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec5->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                        if (part_FTOF_sector_layer2[i] == 6 && part_p[i] > 0 && part_CC_HTCC_nphe[i] > 0)
                            hist_HTCC_Nphe_vs_beta_pip_sec6->Fill(part_CC_HTCC_nphe[i], part_beta[i]);
                    }
                }

                if (FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_DC_z_vertex_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_EC_hit_position_fiducial_check[i])
                {

                    if (part_Cal_PCAL_sector[i] == 1 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec1->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec2->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec3->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec4->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec5->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_Cal_energy_total[i] > 0)
                        hist_sampling_fraction_vs_E_sec6->Fill(part_Cal_energy_total[i], part_Cal_energy_total[i] / part_p[i]);

                    if (part_Cal_PCAL_sector[i] == 1 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec1->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 2 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec2->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 3 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec3->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 4 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec4->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 5 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec5->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                    if (part_Cal_PCAL_sector[i] == 6 && part_p[i] > 0)
                        hist_sampling_fraction_vs_p_sec6->Fill(part_p[i], part_Cal_energy_total[i] / part_p[i]);
                }

            } // end of BUFFER loop for histogram filling

  

            /// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///  fill the statistics:
            /// ///////////////////////////////////////////////////////////////

            hist_electron_count->Fill(e_count);
            hist_proton_count->Fill(p_count);
            hist_neutron_count->Fill(n_count);
            hist_pip_count->Fill(pip_count);
            hist_pim_count->Fill(pim_count);
            hist_Kp_count->Fill(Kp_count);
            hist_Km_count->Fill(Km_count);
            hist_photon_count->Fill(g_count);

        } // end of qa loop

        /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } // end of event loop
    /// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// cut stiatics

    if (show_FD_eid_statistics)
    {
        cout << endl;
        cout << "------------------------------------------------------------------------------------------" << endl;
        cout << "electron PID cut statistics " << endl;
        cout << endl;
        if (FD_eid_default_PID_pass != 0 && FD_eid_charge_pass != 0)
        {
            cout << "number of particles which passed default electron reference PID: " << FD_eid_default_PID_pass << "   ( " << 100 * FD_eid_default_PID_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_default_PID_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed CC nphe cut (not in chain): " << FD_eid_CC_nphe_pass << "   ( " << 100 * FD_eid_CC_nphe_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_CC_nphe_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed EC PCAL vs ECAL cut: " << FD_eid_EC_outer_vs_EC_inner_pass << "   ( " << 100 * FD_eid_EC_outer_vs_EC_inner_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_EC_outer_vs_EC_inner_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed PCAL fiducial cut: " << FD_eid_EC_hit_position_fiducial_pass << "   ( " << 100 * FD_eid_EC_hit_position_fiducial_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_EC_hit_position_fiducial_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed EC sampling fraction cut: " << FD_eid_EC_sampling_fraction_pass << "   ( " << 100 * FD_eid_EC_sampling_fraction_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_EC_sampling_fraction_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_eid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_eid_DC_hit_position_region1_fiducial_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_DC_hit_position_region1_fiducial_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_eid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_eid_DC_hit_position_region2_fiducial_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_DC_hit_position_region2_fiducial_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_eid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_eid_DC_hit_position_region3_fiducial_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_DC_hit_position_region3_fiducial_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed z vertex cut: " << FD_eid_DC_z_vertex_pass << "   ( " << 100 * FD_eid_DC_z_vertex_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_DC_z_vertex_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all FD_eid cuts: " << FD_eid_all_pass << "   ( " << 100 * FD_eid_all_pass / FD_eid_charge_pass << " \% of neg. particles and " << 100 * FD_eid_all_pass / FD_eid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No electrons detected!" << endl;
        cout << endl;
    }

    if (show_charged_ID_statistics)
    {
        cout << "------------------------------------------------------------------------------------------" << endl;
        cout << "charged PID cut statistics " << endl;
        cout << endl;
        cout << "a) Protons " << endl;
        if (FD_protid_default_PID_pass != 0 && FD_protid_charge_pass != 0)
        {
            cout << "number of particles which passed default proton reference PID: " << FD_protid_default_PID_pass << "   ( " << 100 * FD_protid_default_PID_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_default_PID_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_protid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_protid_DC_hit_position_region1_fiducial_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_DC_hit_position_region1_fiducial_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_protid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_protid_DC_hit_position_region2_fiducial_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_DC_hit_position_region2_fiducial_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_protid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_protid_DC_hit_position_region3_fiducial_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_DC_hit_position_region3_fiducial_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_protid_delta_vz_pass << "   ( " << 100 * FD_protid_delta_vz_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_delta_vz_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all proton ID cuts: " << FD_protid_all_pass << "   ( " << 100 * FD_protid_all_pass / FD_protid_charge_pass << " \% of pos. particles and " << 100 * FD_protid_all_pass / FD_protid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No protons detected!" << endl;
        cout << endl;
        cout << "b) Neutron " << endl;
        if (FD_neutrid_default_PID_pass != 0 && FD_neutrid_charge_pass != 0)
        {
            cout << "number of particles which passed default neutron reference PID: " << FD_neutrid_default_PID_pass << "   ( " << 100 * FD_neutrid_default_PID_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_default_PID_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed beta cut: " << FD_neutrid_beta_pass << "   ( " << 100 * FD_neutrid_beta_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_beta_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta beta cut: " << FD_neutrid_delta_beta_pass << "   ( " << 100 * FD_neutrid_delta_beta_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_delta_beta_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed tofmass cut: " << FD_neutrid_tofmass_pass << "   ( " << 100 * FD_neutrid_tofmass_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_tofmass_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_neutrid_delta_vz_pass << "   ( " << 100 * FD_neutrid_delta_vz_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_delta_vz_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all neutron ID cuts: " << FD_neutrid_all_pass << "   ( " << 100 * FD_neutrid_all_pass / FD_neutrid_charge_pass << " \% of neutral particles and " << 100 * FD_neutrid_all_pass / FD_neutrid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No neutrons detected!" << endl;
        cout << endl;
        cout << "b) Pip " << endl;
        if (FD_pipid_default_PID_pass != 0 && FD_pipid_charge_pass != 0)
        {
            cout << "number of particles which passed default pip reference PID: " << FD_pipid_default_PID_pass << "   ( " << 100 * FD_pipid_default_PID_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_default_PID_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_pipid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_pipid_DC_hit_position_region1_fiducial_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_DC_hit_position_region1_fiducial_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_pipid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_pipid_DC_hit_position_region2_fiducial_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_DC_hit_position_region2_fiducial_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_pipid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_pipid_DC_hit_position_region3_fiducial_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_DC_hit_position_region3_fiducial_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_pipid_delta_vz_pass << "   ( " << 100 * FD_pipid_delta_vz_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_delta_vz_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all positive pipon ID cuts: " << FD_pipid_all_pass << "   ( " << 100 * FD_pipid_all_pass / FD_pipid_charge_pass << " \% of pos. particles and " << 100 * FD_pipid_all_pass / FD_pipid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No pi plus detected!" << endl;
        cout << endl;
        cout << "c) Pim " << endl;
        if (FD_pimid_default_PID_pass != 0 && FD_pimid_charge_pass != 0)
        {
            cout << "number of particles which passed default pim reference PID: " << FD_pimid_default_PID_pass << "   ( " << 100 * FD_pimid_default_PID_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_default_PID_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_pimid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_pimid_DC_hit_position_region1_fiducial_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_DC_hit_position_region1_fiducial_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_pimid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_pimid_DC_hit_position_region2_fiducial_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_DC_hit_position_region2_fiducial_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_pimid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_pimid_DC_hit_position_region3_fiducial_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_DC_hit_position_region3_fiducial_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_pimid_delta_vz_pass << "   ( " << 100 * FD_pimid_delta_vz_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_delta_vz_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all negative pion ID cuts: " << FD_pimid_all_pass << "   ( " << 100 * FD_pimid_all_pass / FD_pimid_charge_pass << " \% of neg. particles and " << 100 * FD_pimid_all_pass / FD_pimid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No pi minus detected!" << endl;
        cout << endl;
        cout << "d) Kp " << endl;
        if (FD_Kpid_default_PID_pass != 0 && FD_Kpid_charge_pass != 0)
        {
            cout << "number of particles which passed default Kp reference PID: " << FD_Kpid_default_PID_pass << "   ( " << 100 * FD_Kpid_default_PID_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_default_PID_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_Kpid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_Kpid_DC_hit_position_region1_fiducial_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_DC_hit_position_region1_fiducial_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_Kpid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_Kpid_DC_hit_position_region2_fiducial_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_DC_hit_position_region2_fiducial_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_Kpid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_Kpid_DC_hit_position_region3_fiducial_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_DC_hit_position_region3_fiducial_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_Kpid_delta_vz_pass << "   ( " << 100 * FD_Kpid_delta_vz_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_delta_vz_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all positive Kaon ID cuts: " << FD_Kpid_all_pass << "   ( " << 100 * FD_Kpid_all_pass / FD_Kpid_charge_pass << " \% of pos. particles and " << 100 * FD_Kpid_all_pass / FD_Kpid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No K plus detected!" << endl;
        cout << endl;
        cout << "e) Km " << endl;
        if (FD_Kmid_default_PID_pass != 0 && FD_Kmid_charge_pass != 0)
        {
            cout << "number of particles which passed default Km reference PID: " << FD_Kmid_default_PID_pass << "   ( " << 100 * FD_Kmid_default_PID_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_default_PID_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 1 fid. cut: " << FD_Kmid_DC_hit_position_region1_fiducial_pass << "   ( " << 100 * FD_Kmid_DC_hit_position_region1_fiducial_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_DC_hit_position_region1_fiducial_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 2 fid. cut: " << FD_Kmid_DC_hit_position_region2_fiducial_pass << "   ( " << 100 * FD_Kmid_DC_hit_position_region2_fiducial_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_DC_hit_position_region2_fiducial_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed DC region 3 fid. cut: " << FD_Kmid_DC_hit_position_region3_fiducial_pass << "   ( " << 100 * FD_Kmid_DC_hit_position_region3_fiducial_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_DC_hit_position_region3_fiducial_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed delta vz cut: " << FD_Kmid_delta_vz_pass << "   ( " << 100 * FD_Kmid_delta_vz_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_delta_vz_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all negative Kaon ID cuts: " << FD_Kmid_all_pass << "   ( " << 100 * FD_Kmid_all_pass / FD_Kmid_charge_pass << " \% of neg. particles and " << 100 * FD_Kmid_all_pass / FD_Kmid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No K minus detected!" << endl;
        cout << endl;
    }

    if (show_photon_ID_statistics)
    {
        cout << "------------------------------------------------------------------------------------------" << endl;
        cout << "photon PID cut statistics " << endl;
        cout << endl;
        if (FD_photid_default_PID_pass != 0 && FD_photid_charge_pass != 0)
        {
            cout << "number of particles which passed default photon reference PID: " << FD_photid_default_PID_pass << "   ( " << 100 * FD_photid_default_PID_pass / FD_photid_charge_pass << " \% of neutral particles and " << 100 * FD_photid_default_PID_pass / FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed beta cut: " << FD_photid_beta_pass << "   ( " << 100 * FD_photid_beta_pass / FD_photid_charge_pass << " \% of neutral particles and " << 100 * FD_photid_beta_pass / FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "number of particles which passed EC hit position fiducial cut: " << FD_photid_EC_hit_position_fiducial_pass << "   ( " << 100 * FD_photid_EC_hit_position_fiducial_pass / FD_photid_charge_pass << " \% of neutral particles and " << 100 * FD_photid_EC_hit_position_fiducial_pass / FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
            cout << "__________________________________________________________________________________________________________________________________________________________" << endl;
            cout << "number of particles which passed all photon cuts: " << FD_photid_all_pass << "   ( " << 100 * FD_photid_all_pass / FD_photid_charge_pass << " \% of neutral particles and " << 100 * FD_photid_all_pass / FD_photid_default_PID_pass << " \% of reference PID particles)" << endl;
        }
        else
            cout << "No photons detected!" << endl;
        cout << endl;
    }

    /// /////////////////////////////////////////////////////////////

    cout << endl;
    cout << "Tree successfully analysed!" << endl;
    cout << "Writing the output file ... " << endl;
    out->Write(); // Saving Histograms
    cout << "Histograms saved in File: " << outputfile << endl;
    out->Close(); // Closing Output File
    cout << "... Completed!" << endl;

    return 1;

    /// ///////////////////////////////////////////////////////////////////////////////
} /// end of main
/// ///////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////
/// HTCC

TH2F *create_hist_HTCC_theta_vs_phi(const int cutnum)
{
    char name[100];
    sprintf(name, "HTCC_theta_vs_phi_cut_%02d", cutnum);
    hist_HTCC_theta_vs_phi[cutnum] = new TH2F(name, name, 360, -180, 180, 80, 0, 40);
    hist_HTCC_theta_vs_phi[cutnum]->GetXaxis()->SetTitle("#Phi_{CC}");
    hist_HTCC_theta_vs_phi[cutnum]->GetYaxis()->SetTitle("#theta_{CC}");
    return (hist_HTCC_theta_vs_phi[cutnum]);
}

TH1F *create_hist_HTCC_nphe(const int cutnum)
{
    char name[100];
    sprintf(name, "HTCC_nphe_cut_%02d", cutnum);
    hist_HTCC_nphe[cutnum] = new TH1F(name, name, 70, 0, 70);
    hist_HTCC_nphe[cutnum]->GetXaxis()->SetTitle("nphe");
    hist_HTCC_nphe[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_HTCC_nphe[cutnum]);
}

TH2F *create_hist_HTCC_nphe_vs_sampling_fraction(const int cutnum)
{
    char name[100];
    sprintf(name, "HTCC_nphe_vs_sampfrac_cut_%02d", cutnum);
    hist_HTCC_nphe_vs_sampling_fraction[cutnum] = new TH2F(name, name, 70, 0, 70, 60, 0, 0.6);
    hist_HTCC_nphe_vs_sampling_fraction[cutnum]->GetXaxis()->SetTitle("nphe");
    hist_HTCC_nphe_vs_sampling_fraction[cutnum]->GetYaxis()->SetTitle("EC sampling fraction");
    return (hist_HTCC_nphe_vs_sampling_fraction[cutnum]);
}

// ////////////////////////////////////////////////////////////////////////////////////
// EC cuts

TH2F *create_hist_EC_PCAL_vs_EC_ECAL(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_vs_EC_ECAL_cut_%02d", cutnum);
    hist_EC_PCAL_vs_EC_ECAL[cutnum] = new TH2F(name, name, 1000, 0, 1.0, 1000, 0, 1.0);
    hist_EC_PCAL_vs_EC_ECAL[cutnum]->GetXaxis()->SetTitle("E_{PCAL} /GeV");
    hist_EC_PCAL_vs_EC_ECAL[cutnum]->GetYaxis()->SetTitle("E_{ECAL} /GeV");
    return (hist_EC_PCAL_vs_EC_ECAL[cutnum]);
}

TH2F *create_hist_EC_outer_vs_EC_inner(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_outer_vs_EC_inner_cut_%02d", cutnum);
    hist_EC_outer_vs_EC_inner[cutnum] = new TH2F(name, name, 1000, 0, 1.0, 1000, 0, 1.0);
    hist_EC_outer_vs_EC_inner[cutnum]->GetXaxis()->SetTitle("E_{EC inner} /GeV");
    hist_EC_outer_vs_EC_inner[cutnum]->GetYaxis()->SetTitle("E_{EC outer} /GeV");
    return (hist_EC_outer_vs_EC_inner[cutnum]);
}

///

TH2F *create_hist_EC_total_sampling_fraction_sec1(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec1_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec1[cutnum]);
}
TH2F *create_hist_EC_total_sampling_fraction_sec2(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec2_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec2[cutnum]);
}
TH2F *create_hist_EC_total_sampling_fraction_sec3(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec3_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec3[cutnum]);
}
TH2F *create_hist_EC_total_sampling_fraction_sec4(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec4_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec4[cutnum]);
}
TH2F *create_hist_EC_total_sampling_fraction_sec5(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec5_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec5[cutnum]);
}
TH2F *create_hist_EC_total_sampling_fraction_sec6(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_total_sampling_fraction_sec6_cut_%02d", cutnum);
    hist_EC_total_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100, 0, 11, 600, 0, 0.6);
    hist_EC_total_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_total_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_total_sampling_fraction_sec6[cutnum]);
}

TH2F *create_hist_EC_PCAL_sampling_fraction_sec1(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec1_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec1[cutnum]);
}
TH2F *create_hist_EC_PCAL_sampling_fraction_sec2(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec2_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec2[cutnum]);
}
TH2F *create_hist_EC_PCAL_sampling_fraction_sec3(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec3_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec3[cutnum]);
}
TH2F *create_hist_EC_PCAL_sampling_fraction_sec4(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec4_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec4[cutnum]);
}
TH2F *create_hist_EC_PCAL_sampling_fraction_sec5(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec5_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec5[cutnum]);
}
TH2F *create_hist_EC_PCAL_sampling_fraction_sec6(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_sampling_fraction_sec6_cut_%02d", cutnum);
    hist_EC_PCAL_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100, 0, 11, 60, 0, 0.6);
    hist_EC_PCAL_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_PCAL_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_PCAL_sampling_fraction_sec6[cutnum]);
}

TH2F *create_hist_EC_ECAL_sampling_fraction_sec1(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec1_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec1[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec1[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec1[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec1[cutnum]);
}
TH2F *create_hist_EC_ECAL_sampling_fraction_sec2(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec2_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec2[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec2[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec2[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec2[cutnum]);
}
TH2F *create_hist_EC_ECAL_sampling_fraction_sec3(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec3_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec3[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec3[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec3[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec3[cutnum]);
}
TH2F *create_hist_EC_ECAL_sampling_fraction_sec4(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec4_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec4[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec4[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec4[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec4[cutnum]);
}
TH2F *create_hist_EC_ECAL_sampling_fraction_sec5(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec5_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec5[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec5[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec5[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec5[cutnum]);
}
TH2F *create_hist_EC_ECAL_sampling_fraction_sec6(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_ECAL_sampling_fraction_sec6_cut_%02d", cutnum);
    hist_EC_ECAL_sampling_fraction_sec6[cutnum] = new TH2F(name, name, 1100, 0, 11, 40, 0, 0.4);
    hist_EC_ECAL_sampling_fraction_sec6[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_ECAL_sampling_fraction_sec6[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_ECAL_sampling_fraction_sec6[cutnum]);
}

///

TH2F *create_hist_EC_PCAL_hit_position(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_hit_position_cut_%02d", cutnum);
    hist_EC_PCAL_hit_position[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_EC_PCAL_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_EC_PCAL_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_EC_PCAL_hit_position[cutnum]);
}

TH2F *create_hist_EC_inner_hit_position(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_inner_hit_position_cut_%02d", cutnum);
    hist_EC_inner_hit_position[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_EC_inner_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_EC_inner_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_EC_inner_hit_position[cutnum]);
}

TH2F *create_hist_EC_outer_hit_position(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_outer_hit_position_cut_%02d", cutnum);
    hist_EC_outer_hit_position[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_EC_outer_hit_position[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_EC_outer_hit_position[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_EC_outer_hit_position[cutnum]);
}

// ////////////////////////////////////////////////////////////////////////////////////////
// DC cuts

TH2F *create_hist_DC_hit_position_region1(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region1_cut_%02d", cutnum);
    hist_DC_hit_position_region1[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_DC_hit_position_region1[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region1[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region1[cutnum]);
}

TH2F *create_hist_DC_hit_position_region2(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region2_cut_%02d", cutnum);
    hist_DC_hit_position_region2[cutnum] = new TH2F(name, name, 1000, -500, 500, 1000, -500, 500);
    hist_DC_hit_position_region2[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region2[cutnum]);
}

TH2F *create_hist_DC_hit_position_region3(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region3_cut_%02d", cutnum);
    hist_DC_hit_position_region3[cutnum] = new TH2F(name, name, 1000, -500, 500, 1000, -500, 500);
    hist_DC_hit_position_region3[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region3[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region3[cutnum]);
}

TH1F *create_hist_DC_z_vertex_sec1(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec1_cut_%02d", cutnum);
    hist_DC_z_vertex_sec1[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec1[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec1[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec1[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec2(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec2_cut_%02d", cutnum);
    hist_DC_z_vertex_sec2[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec2[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec2[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec2[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec3(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec3_cut_%02d", cutnum);
    hist_DC_z_vertex_sec3[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec3[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec3[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec3[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec4(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec4_cut_%02d", cutnum);
    hist_DC_z_vertex_sec4[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec4[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec4[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec4[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec5(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec5_cut_%02d", cutnum);
    hist_DC_z_vertex_sec5[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec5[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec5[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec5[cutnum]);
}
TH1F *create_hist_DC_z_vertex_sec6(int cutnum)
{
    char name[100];
    sprintf(name, "DC_z_vertex_sec6_cut_%02d", cutnum);
    hist_DC_z_vertex_sec6[cutnum] = new TH1F(name, name, 400, -40, 40);
    hist_DC_z_vertex_sec6[cutnum]->GetXaxis()->SetTitle("v_{z} /cm");
    hist_DC_z_vertex_sec6[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_DC_z_vertex_sec6[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// Basic hadron cuts

TH2F *create_DC_hit_position_region1_hadron(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region1_hadron_cut_%02d", cutnum);
    hist_DC_hit_position_region1_hadron[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_DC_hit_position_region1_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region1_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region1_hadron[cutnum]);
}

TH2F *create_DC_hit_position_region2_hadron(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region2_hadron_cut_%02d", cutnum);
    hist_DC_hit_position_region2_hadron[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_DC_hit_position_region2_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region2_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region2_hadron[cutnum]);
}

TH2F *create_DC_hit_position_region3_hadron(int cutnum)
{
    char name[100];
    sprintf(name, "DC_hit_position_region3_hadron_cut_%02d", cutnum);
    hist_DC_hit_position_region3_hadron[cutnum] = new TH2F(name, name, 900, -450, 450, 900, -450, 450);
    hist_DC_hit_position_region3_hadron[cutnum]->GetXaxis()->SetTitle("x /cm");
    hist_DC_hit_position_region3_hadron[cutnum]->GetYaxis()->SetTitle("y /cm");
    return (hist_DC_hit_position_region3_hadron[cutnum]);
}

TH2F *create_hist_EC_outer_vs_EC_inner_hadron(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_outer_vs_EC_inner_hadron_cut_%02d", cutnum);
    hist_EC_outer_vs_EC_inner_hadron[cutnum] = new TH2F(name, name, 1000, 0, 1.0, 1000, 0, 1.0);
    hist_EC_outer_vs_EC_inner_hadron[cutnum]->GetXaxis()->SetTitle("E_{EC inner} /GeV");
    hist_EC_outer_vs_EC_inner_hadron[cutnum]->GetYaxis()->SetTitle("E_{EC outer} /GeV");
    return (hist_EC_outer_vs_EC_inner_hadron[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// TOF beta cuts

TH2F *create_hist_beta_vs_p(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_cut_%02d", cutnum);
    hist_beta_vs_p[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p[cutnum]);
}

TH2F *create_hist_beta_vs_p_sec1(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec1_cut_%02d", cutnum);
    hist_beta_vs_p_sec1[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec1[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec1[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec1[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec2(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec2_cut_%02d", cutnum);
    hist_beta_vs_p_sec2[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec2[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec2[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec2[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec3(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec3_cut_%02d", cutnum);
    hist_beta_vs_p_sec3[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec3[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec3[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec3[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec4(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec4_cut_%02d", cutnum);
    hist_beta_vs_p_sec4[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec4[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec4[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec4[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec5(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec5_cut_%02d", cutnum);
    hist_beta_vs_p_sec5[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec5[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec5[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec5[cutnum]);
}
TH2F *create_hist_beta_vs_p_sec6(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_momentum_sec6_cut_%02d", cutnum);
    hist_beta_vs_p_sec6[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_beta_vs_p_sec6[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_beta_vs_p_sec6[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_sec6[cutnum]);
}

TH2F *create_hist_delta_beta_vs_p(int cutnum)
{
    char name[100];
    sprintf(name, "delta_beta_vs_momentum_cut_%02d", cutnum);
    hist_delta_beta_vs_p[cutnum] = new TH2F(name, name, 1100, 0, 11, 250, -0.5, 0.5);
    hist_delta_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_delta_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#Delta #beta");
    return (hist_delta_beta_vs_p[cutnum]);
}

TH1F *create_hist_delta_beta(int cutnum)
{
    char name[100];
    sprintf(name, "delta_beta_cut_%02d", cutnum);
    hist_delta_beta[cutnum] = new TH1F(name, name, 500, -0.5, 0.5);
    hist_delta_beta[cutnum]->GetXaxis()->SetTitle("#Delta #beta");
    hist_delta_beta[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_delta_beta[cutnum]);
}

TH2F *create_hist_tofmass_vs_p(int cutnum)
{
    char name[100];
    sprintf(name, "tofmass2_vs_momentum_cut_%02d", cutnum);
    hist_tofmass_vs_p[cutnum] = new TH2F(name, name, 1100, 0, 11, 1100, -0.2, 2.0);
    hist_tofmass_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_tofmass_vs_p[cutnum]->GetYaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}");
    return (hist_tofmass_vs_p[cutnum]);
}

TH1F *create_hist_tofmass(int cutnum)
{
    char name[100];
    sprintf(name, "tofmass2_cut_%02d", cutnum);
    hist_tofmass[cutnum] = new TH1F(name, name, 1100, -0.2, 2.0);
    hist_tofmass[cutnum]->GetXaxis()->SetTitle("m_{TOF}^{2} /GeV^{2}");
    hist_tofmass[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_tofmass[cutnum]);
}

TH1F *create_hist_delta_vz(int cutnum)
{
    char name[100];
    sprintf(name, "delta_vz_cut_%02d", cutnum);
    hist_delta_vz[cutnum] = new TH1F(name, name, 200, -100, 100);
    hist_delta_vz[cutnum]->GetXaxis()->SetTitle("#Delta v_{z} /cm");
    hist_delta_vz[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_delta_vz[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// CTOF beta plots

TH2F *create_hist_CD_beta_vs_p(int cutnum)
{
    char name[100];
    sprintf(name, "CD_beta_vs_momentum_cut_%02d", cutnum);
    hist_CD_beta_vs_p[cutnum] = new TH2F(name, name, 1100, 0, 11, 1400, 0, 1.40);
    hist_CD_beta_vs_p[cutnum]->GetXaxis()->SetTitle("p /GeV");
    hist_CD_beta_vs_p[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_CD_beta_vs_p[cutnum]);
}

TH1F *create_hist_CD_delta_vz(int cutnum)
{
    char name[100];
    sprintf(name, "CD_delta_vz_cut_%02d", cutnum);
    hist_CD_delta_vz[cutnum] = new TH1F(name, name, 200, -100, 100);
    hist_CD_delta_vz[cutnum]->GetXaxis()->SetTitle("#Delta v_{z} /cm");
    hist_CD_delta_vz[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_CD_delta_vz[cutnum]);
}

// /////////////////////////////////////////////////////////////////////////////
// photon ID plots based on ECAL

TH2F *create_hist_beta_vs_p_phot(int cutnum)
{
    char name[100];
    sprintf(name, "beta_vs_p_phot_cut_%02d", cutnum);
    hist_beta_vs_p_phot[cutnum] = new TH2F(name, name, 300, 0, 6, 700, 0, 1.40);
    hist_beta_vs_p_phot[cutnum]->GetXaxis()->SetTitle("p / GeV");
    hist_beta_vs_p_phot[cutnum]->GetYaxis()->SetTitle("#beta");
    return (hist_beta_vs_p_phot[cutnum]);
}

TH1F *create_hist_beta_phot(int cutnum)
{
    char name[100];
    sprintf(name, "beta_phot_cut_%02d", cutnum);
    hist_beta_phot[cutnum] = new TH1F(name, name, 750, 0, 1.5);
    hist_beta_phot[cutnum]->GetXaxis()->SetTitle("#beta");
    hist_beta_phot[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_beta_phot[cutnum]);
}

TH2F *create_hist_EC_sampling_fraction_phot(int cutnum)
{
    char name[100];
    sprintf(name, "EC_sampling_fraction_phot_cut_%02d", cutnum);
    hist_EC_sampling_fraction_phot[cutnum] = new TH2F(name, name, 600, 0, 6, 100, 0, 1);
    hist_EC_sampling_fraction_phot[cutnum]->GetXaxis()->SetTitle("momentum");
    hist_EC_sampling_fraction_phot[cutnum]->GetYaxis()->SetTitle("sampling fraction");
    return (hist_EC_sampling_fraction_phot[cutnum]);
}

TH2F *create_hist_EC_PCAL_vs_EC_ECAL_phot(const int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_vs_EC_ECAL_phot_cut_%02d", cutnum);
    hist_EC_PCAL_vs_EC_ECAL_phot[cutnum] = new TH2F(name, name, 400, 0, 0.4, 400, 0, 0.4);
    hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]->GetXaxis()->SetTitle("E_{PCAL} /GeV");
    hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]->GetYaxis()->SetTitle("E_{ECAL} /GeV");
    return (hist_EC_PCAL_vs_EC_ECAL_phot[cutnum]);
}

TH2F *create_hist_EC_PCAL_hit_position_phot(int cutnum)
{
    char name[100];
    sprintf(name, "EC_PCAL_hit_position_phot_cut_%02d", cutnum);
    hist_EC_PCAL_hit_position_phot[cutnum] = new TH2F(name, name, 500, -500, 500, 500, -500, 500);
    hist_EC_PCAL_hit_position_phot[cutnum]->GetXaxis()->SetTitle("PCAL x /cm");
    hist_EC_PCAL_hit_position_phot[cutnum]->GetYaxis()->SetTitle("PCAL y /cm");
    return (hist_EC_PCAL_hit_position_phot[cutnum]);
}

/// ////////////////////////////////////////////////////////////////////////
/// FT PID plots

TH2F *create_hist_FT_FTCAL_energy_vs_radius(int cutnum)
{
    char name[100];
    sprintf(name, "FT_FTCAL_energy_vs_radius_cut_%02d", cutnum);
    hist_FT_FTCAL_energy_vs_radius[cutnum] = new TH2F(name, name, 600, 0, 6, 100, 0, 10);
    hist_FT_FTCAL_energy_vs_radius[cutnum]->GetXaxis()->SetTitle("FTCAL energy /GeV");
    hist_FT_FTCAL_energy_vs_radius[cutnum]->GetYaxis()->SetTitle("FTCAL cluster radius /cm");
    return (hist_FT_FTCAL_energy_vs_radius[cutnum]);
}

TH2F *create_hist_FT_FTCAL_hit_position(int cutnum)
{
    char name[100];
    sprintf(name, "FT_FTCAL_hit_position_cut_%02d", cutnum);
    hist_FT_FTCAL_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);
    hist_FT_FTCAL_hit_position[cutnum]->GetXaxis()->SetTitle("FTCAL x /cm");
    hist_FT_FTCAL_hit_position[cutnum]->GetYaxis()->SetTitle("FTCAL y /cm");
    return (hist_FT_FTCAL_hit_position[cutnum]);
}

TH2F *create_hist_FT_FTTRK_hit_position(int cutnum)
{
    char name[100];
    sprintf(name, "FT_FTTRK_hit_position_cut_%02d", cutnum);
    hist_FT_FTTRK_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);
    hist_FT_FTTRK_hit_position[cutnum]->GetXaxis()->SetTitle("FTTRK x /cm");
    hist_FT_FTTRK_hit_position[cutnum]->GetYaxis()->SetTitle("FTTRK y /cm");
    return (hist_FT_FTTRK_hit_position[cutnum]);
}

TH2F *create_hist_FT_FTHODO_hit_position(int cutnum)
{
    char name[100];
    sprintf(name, "FT_FTHODO_hit_position_cut_%02d", cutnum);
    hist_FT_FTHODO_hit_position[cutnum] = new TH2F(name, name, 40, -20, 20, 40, -20, 20);
    hist_FT_FTHODO_hit_position[cutnum]->GetXaxis()->SetTitle("FTHODO x /cm");
    hist_FT_FTHODO_hit_position[cutnum]->GetYaxis()->SetTitle("FTHODO y /cm");
    return (hist_FT_FTHODO_hit_position[cutnum]);
}

TH1F *create_hist_FT_beta(int cutnum)
{
    char name[100];
    sprintf(name, "FT_beta_cut_%02d", cutnum);
    hist_FT_beta[cutnum] = new TH1F(name, name, 130, 0, 1.3);
    hist_FT_beta[cutnum]->GetXaxis()->SetTitle("#beta");
    hist_FT_beta[cutnum]->GetYaxis()->SetTitle("counts");
    return (hist_FT_beta[cutnum]);
}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////
/// assign raw particles

void get_event_properties(event_ptr event, mcevt_ptr mcevent)
{
    if(simulation)
    {
        //MC_helicity = mcevent->getGenHelicity(); //not available
        MC_Ebeam = mcevent->getEbeam();
        MC_weight = mcevent->getWeight();
        MC_Npart = mcevent->getNpart();
    }

    STTime = event->getStartTime();
    RFTime = event->getRFTime();
    Helic = event->getHelicity();
    Helic_raw = event->getHelicityRaw();
    beam_charge = event->getBeamCharge();
}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////
/// assign raw particles and detector properties

void assign_particles(std::vector<region_part_ptr>& particles, mcpar_ptr mcparticles)
{
    Npart = particles.size();

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {
            part_px[i] = particles[i]->par()->getPx();
            part_py[i] = particles[i]->par()->getPy();
            part_pz[i] = particles[i]->par()->getPz();
            part_vx[i] = particles[i]->par()->getVx();
            part_vy[i] = particles[i]->par()->getVy();
            part_vz[i] = particles[i]->par()->getVz();
            part_charge[i] = particles[i]->par()->getCharge();
            part_beta[i] = particles[i]->par()->getBeta();
            part_pid[i] = particles[i]->par()->getPid();
            part_status[i] = particles[i]->par()->getStatus();
            part_chi2pid[i] = particles[i]->par()->getChi2Pid();
            part_p[i] = particles[i]->par()->getP();
            part_theta[i] = particles[i]->getTheta();
            part_phi[i] = particles[i]->getPhi();
        }
    }


    int Cal_Nentries = 0;
    int CC_Nentries = 0;
    int FT_Nentries = 0;
    int SC_Nentries = 0;
    int Traj_Nentries = 0;
    int TRK_Nentries = 0;
    int TBT_Nentries = 0;
    int RICHHadCher_Nentries = 0;
    int RICHHadrons_Nentries = 0;

    Cal_Nentries = Npart;
    CC_Nentries = Npart;
    FT_Nentries = Npart;
    SC_Nentries = Npart;
    TRK_Nentries = Npart;
    Traj_Nentries = Npart;
    if (userich)
        RICHHadCher_Nentries = Npart;
    if (userich)
        RICHHadrons_Nentries = Npart;

    // Calorimeter bank (detector = 7  ---  layer:  PCAL = 1, ECin = 4, ECout = 7)
    

    if (Cal_Nentries > 0)
    {
        for (int i = 0; i < Cal_Nentries; i++)
        {
            part_Cal_PCAL_sector[i] = particles[i]->cal(PCAL)->getSector();
            part_Cal_PCAL_energy[i] = particles[i]->cal(PCAL)->getEnergy();
            part_Cal_PCAL_time[i] = particles[i]->cal(PCAL)->getTime();
            part_Cal_PCAL_path[i] = particles[i]->cal(PCAL)->getPath();
            part_Cal_PCAL_x[i] = particles[i]->cal(PCAL)->getX();
            part_Cal_PCAL_y[i] = particles[i]->cal(PCAL)->getY();
            part_Cal_PCAL_z[i] = particles[i]->cal(PCAL)->getZ();
            part_Cal_PCAL_lu[i] = particles[i]->cal(PCAL)->getLu();
            part_Cal_PCAL_lv[i] = particles[i]->cal(PCAL)->getLv();
            part_Cal_PCAL_lw[i] = particles[i]->cal(PCAL)->getLw();

            part_Cal_PCAL_m2u[i] = particles[i]->cal(PCAL)->getM2u();
            part_Cal_PCAL_m2v[i] = particles[i]->cal(PCAL)->getM2v();
            part_Cal_PCAL_m2w[i] = particles[i]->cal(PCAL)->getM2w();
            part_Cal_PCAL_m3u[i] = particles[i]->cal(PCAL)->getM3u();
            part_Cal_PCAL_m3v[i] = particles[i]->cal(PCAL)->getM3v();
            part_Cal_PCAL_m3w[i] = particles[i]->cal(PCAL)->getM3w();

            part_Cal_ECin_sector[i] = particles[i]->cal(ECIN)->getSector();
            part_Cal_ECin_energy[i] = particles[i]->cal(ECIN)->getEnergy();
            part_Cal_ECin_time[i] = particles[i]->cal(ECIN)->getTime();
            part_Cal_ECin_path[i] = particles[i]->cal(ECIN)->getPath();
            part_Cal_ECin_x[i] = particles[i]->cal(ECIN)->getX();
            part_Cal_ECin_y[i] = particles[i]->cal(ECIN)->getY();
            part_Cal_ECin_z[i] = particles[i]->cal(ECIN)->getZ();
            part_Cal_ECin_lu[i] = particles[i]->cal(ECIN)->getLu();
            part_Cal_ECin_lv[i] = particles[i]->cal(ECIN)->getLv();
            part_Cal_ECin_lw[i] = particles[i]->cal(ECIN)->getLw();

            part_Cal_ECin_m2u[i] = particles[i]->cal(ECIN)->getM2u();
            part_Cal_ECin_m2v[i] = particles[i]->cal(ECIN)->getM2v();
            part_Cal_ECin_m2w[i] = particles[i]->cal(ECIN)->getM2w();
            part_Cal_ECin_m3u[i] = particles[i]->cal(ECIN)->getM3u();
            part_Cal_ECin_m3v[i] = particles[i]->cal(ECIN)->getM3v();
            part_Cal_ECin_m3w[i] = particles[i]->cal(ECIN)->getM3w();

            part_Cal_ECout_sector[i] = particles[i]->cal(ECOUT)->getSector();
            part_Cal_ECout_energy[i] = particles[i]->cal(ECOUT)->getEnergy();
            part_Cal_ECout_time[i] = particles[i]->cal(ECOUT)->getTime();
            part_Cal_ECout_path[i] = particles[i]->cal(ECOUT)->getPath();
            part_Cal_ECout_x[i] = particles[i]->cal(ECOUT)->getX();
            part_Cal_ECout_y[i] = particles[i]->cal(ECOUT)->getY();
            part_Cal_ECout_z[i] = particles[i]->cal(ECOUT)->getZ();
            part_Cal_ECout_lu[i] = particles[i]->cal(ECOUT)->getLu();
            part_Cal_ECout_lv[i] = particles[i]->cal(ECOUT)->getLv();
            part_Cal_ECout_lw[i] = particles[i]->cal(ECOUT)->getLw();

            part_Cal_ECout_m2u[i] = particles[i]->cal(ECOUT)->getM2u();
            part_Cal_ECout_m2v[i] = particles[i]->cal(ECOUT)->getM2v();
            part_Cal_ECout_m2w[i] = particles[i]->cal(ECOUT)->getM2w();
            part_Cal_ECout_m3u[i] = particles[i]->cal(ECOUT)->getM3u();
            part_Cal_ECout_m3v[i] = particles[i]->cal(ECOUT)->getM3v();
            part_Cal_ECout_m3w[i] = particles[i]->cal(ECOUT)->getM3w();
        }
    }


    for (Int_t i = 0; i < BUFFER; i++)
    {
        part_Cal_energy_total[i] = part_Cal_PCAL_energy[i] + part_Cal_ECin_energy[i] + part_Cal_ECout_energy[i];
    }

    // Cherenkov bank (detectors:  HTCC = 15,  LTCC = 16)

    if (CC_Nentries > 0)
    {
        for (int i = 0; i < CC_Nentries; i++)
        {
            part_CC_HTCC_sector[i] = particles[i]->che(HTCC)->getSector();
            part_CC_HTCC_nphe[i] = particles[i]->che(HTCC)->getNphe();
            part_CC_HTCC_time[i] = particles[i]->che(HTCC)->getTime();
            part_CC_HTCC_path[i] = particles[i]->che(HTCC)->getPath();
            part_CC_HTCC_theta[i] = particles[i]->che(HTCC)->getDtheta();
            part_CC_HTCC_phi[i] = particles[i]->che(HTCC)->getDPhi();

            part_CC_LTCC_sector[i] = particles[i]->che(LTCC)->getSector();
            part_CC_LTCC_nphe[i] = particles[i]->che(LTCC)->getNphe();
            part_CC_LTCC_time[i] = particles[i]->che(LTCC)->getTime();
            part_CC_LTCC_path[i] = particles[i]->che(LTCC)->getPath();
            part_CC_LTCC_theta[i] = particles[i]->che(LTCC)->getDtheta();
            part_CC_LTCC_phi[i] = particles[i]->che(LTCC)->getDPhi();
        }
    }

    // Forward Tagger bank (detectors: FTCAL = 10, FTTRK = 13, FTHODO = 11)

    if (FT_Nentries > 0)
    {
        for (int i = 0; i < FT_Nentries; i++)
        {
            part_FT_energy[i] = particles[i]->ft(FTCAL)->getEnergy();
            part_FT_radius[i] = particles[i]->ft(FTCAL)->getRadius();
            part_FTHODO_time[i] = particles[i]->ft(FTHODO)->getTime();
            part_FTHODO_path[i] = particles[i]->ft(FTHODO)->getPath();
            part_FTCAL_time[i] = particles[i]->ft(FTCAL)->getTime();
            part_FTCAL_path[i] = particles[i]->ft(FTCAL)->getPath();
            part_FTCAL_x[i] = particles[i]->ft(FTCAL)->getX();
            part_FTCAL_y[i] = particles[i]->ft(FTCAL)->getY();
            part_FTCAL_z[i] = particles[i]->ft(FTCAL)->getZ();
            part_FTHODO_x[i] = particles[i]->ft(FTHODO)->getX();
            part_FTHODO_y[i] = particles[i]->ft(FTHODO)->getY();
            part_FTHODO_z[i] = particles[i]->ft(FTHODO)->getZ();
            part_FTTRK_x[i] = particles[i]->ft(FTTRK)->getX(); // DET LAYER ID 13 - FTTRK
            part_FTTRK_y[i] = particles[i]->ft(FTTRK)->getY();
            part_FTTRK_z[i] = particles[i]->ft(FTTRK)->getZ();
        }
    }

    // Scintillator bank (detectors:  CTOF = 4,  CND = 3,  FTOF = 12  ---  for FTOF:  layer = 1,2)

    if (SC_Nentries > 0)
    {
        for (int i = 0; i < SC_Nentries; i++)
        {
            // FTOF

            part_FTOF_layer[i] = particles[i]->sci(FTOF1A)->getLayer();
            part_FTOF_energy_layer1[i] = particles[i]->sci(FTOF1A)->getEnergy();
            part_FTOF_time_layer1[i] = particles[i]->sci(FTOF1A)->getTime();
            part_FTOF_path_layer1[i] = particles[i]->sci(FTOF1A)->getPath();
            part_FTOF_sector_layer1[i] = particles[i]->sci(FTOF1A)->getSector();
            part_FTOF_component_layer1[i] = particles[i]->sci(FTOF1A)->getComponent();

            part_FTOF_layer[i] = particles[i]->sci(FTOF1B)->getLayer();
            part_FTOF_sector_layer2[i] = particles[i]->sci(FTOF1B)->getSector();
            part_FTOF_component_layer2[i] = particles[i]->sci(FTOF1B)->getComponent();

            part_FTOF_layer[i] = particles[i]->sci(FTOF2)->getLayer();
            part_FTOF_energy_layer3[i] = particles[i]->sci(FTOF2)->getEnergy();
            part_FTOF_time_layer3[i] = particles[i]->sci(FTOF2)->getTime();
            part_FTOF_path_layer3[i] = particles[i]->sci(FTOF2)->getPath();
            part_FTOF_sector_layer3[i] = particles[i]->sci(FTOF2)->getSector();
            part_FTOF_component_layer3[i] = particles[i]->sci(FTOF2)->getComponent();

            part_FTOF_energy[i] = particles[i]->sci(FTOF1B)->getEnergy();
            part_FTOF_time[i] = particles[i]->sci(FTOF1B)->getTime();
            part_FTOF_path[i] = particles[i]->sci(FTOF1B)->getPath();

            part_CTOF_component[i] = particles[i]->sci(CTOF)->getComponent();
            part_CTOF_energy[i] = particles[i]->sci(CTOF)->getEnergy();
            part_CTOF_time[i] = particles[i]->sci(CTOF)->getTime();
            part_CTOF_path[i] = particles[i]->sci(CTOF)->getPath();

            part_CND_component[i] = particles[i]->sci(CND2)->getComponent();
            part_CND_energy[i] = particles[i]->sci(CND2)->getEnergy();
            part_CND_time[i] = particles[i]->sci(CND2)->getTime();
            part_CND_path[i] = particles[i]->sci(CND2)->getPath();
        }
    }

    // tracking banks  (detectors: DC = 6, BST = 2,  BMT = 1, FMT = 8)

    if (TRK_Nentries > 0)
    {
        for (int i = 0; i < TRK_Nentries; i++)
        {
                part_DC_Track_chi2[i] = particles[i]->trk(DC)->getChi2();
                part_DC_Track_NDF[i] = particles[i]->trk(DC)->getNDF();
                part_DC_Track_status[i] = particles[i]->trk(DC)->getStatus();
        }
    }

    // trajectory crosses  (layers: 6 = DC region 1 start,  18 = DC region 2 center,  36 = DC region 3 end )

    if (Traj_Nentries > 0)
    {
        for (int i = 0; i < Traj_Nentries; i++)
        {
                part_DC_c1x[i] = particles[i]->traj(DC, 6)->getX();
                part_DC_c1y[i] = particles[i]->traj(DC, 6)->getY();
                part_DC_c1z[i] = particles[i]->traj(DC, 6)->getZ();

                part_DC_c2x[i] = particles[i]->traj(DC, 18)->getX();
                part_DC_c2y[i] = particles[i]->traj(DC, 18)->getY();
                part_DC_c2z[i] = particles[i]->traj(DC, 18)->getZ();

                part_DC_c3x[i] = particles[i]->traj(DC, 36)->getX();
                part_DC_c3y[i] = particles[i]->traj(DC, 36)->getY();
                part_DC_c3z[i] = particles[i]->traj(DC, 36)->getZ();
                
                part_DC_edge1[i]= particles[i]->traj(DC, 6)->getEdge();
                part_DC_edge2[i]= particles[i]->traj(DC, 18)->getEdge();
                part_DC_edge3[i]= particles[i]->traj(DC, 36)->getEdge();
        }
    }

    for (int i = 0; i < BUFFER; ++i)
    {
        part_DC_sector[i] = determineSector(i);
    }

    if (userich)
    {
        if (RICHHadCher_Nentries > 0)
        {
            for (int i = 0; i < RICHHadCher_Nentries; i++)
            {
                    part_RICH_best[i] = particles[i]->rich()->getBest_PID();
            }
        }
    }
    
    /// //////////////////////////////////////////////
    /// MC:

    if (simulation == true && mcparticles != nullptr)
    {

        int MC_Particle_Nentries = 0;
        int MC_Lund_Nentries = 0;

        MC_Particle_Nentries = MC_Npart;
        MC_Lund_Nentries = MC_Npart;


        if (MC_Particle_Nentries > 0)
        {
            for (int i = 0; i < MC_Particle_Nentries; i++)
            {
                if (i < BUFFER)
                {
                    partMC_pid[i] = mcparticles->getPid(i);
                    partMC_px[i] = mcparticles->getPx(i);
                    partMC_py[i] = mcparticles->getPy(i);
                    partMC_pz[i] = mcparticles->getPz(i);
                    partMC_vx[i] = mcparticles->getVx(i);
                    partMC_vy[i] = mcparticles->getVy(i);
                    partMC_vz[i] = mcparticles->getVz(i);
                    partMC_mother[i] = mcparticles->getParent(i);
                    partMC_p[i] = sqrt(partMC_px[i] * partMC_px[i] + partMC_py[i] * partMC_py[i] + partMC_pz[i] * partMC_pz[i]);
                    partMC_theta[i] = acos(partMC_pz[i] / partMC_p[i]);
                    partMC_phi[i] = atan2(partMC_py[i], partMC_px[i]);
                    
                     
                    /*
                    partMC_pid[i] = particles[i]->mc()->getPid();
                    partMC_px[i] = particles[i]->mc()->getPx();
                    partMC_py[i] = particles[i]->mc()->getPy();
                    partMC_pz[i] = particles[i]->mc()->getPz();
                    partMC_vx[i] = particles[i]->mc()->getVx();
                    partMC_vy[i] = particles[i]->mc()->getVy();
                    partMC_vz[i] = particles[i]->mc()->getVz();
                    partMC_p[i] = sqrt(partMC_px[i] * partMC_px[i] + partMC_py[i] * partMC_py[i] + partMC_pz[i] * partMC_pz[i]);
                    partMC_theta[i] = acos(partMC_pz[i] / partMC_p[i]);
                    partMC_phi[i] = atan2(partMC_py[i], partMC_px[i]);
                    */
                
                }
            }
        }
    
        /*
        if (MC_Lund_Nentries > 0)
        {
            for (int i = 0; i < MC_Lund_Nentries; i++)
            {
                if (i < BUFFER)
                {
                    partLUND_pid[i] = mcparticles->getPid(i);
                    partLUND_mass[i] = mcparticles->getMass(i);
                    partLUND_px[i] = mcparticles->getPx(i);
                    partLUND_py[i] = mcparticles->getPy(i);
                    partLUND_pz[i] = mcparticles->getPz(i);
                    partLUND_vx[i] = mcparticles->getVx(i);
                    partLUND_vy[i] = mcparticles->getVy(i);
                    partLUND_vz[i] = mcparticles->getVz(i);
                    partLUND_mother[i] = mcparticles->getParent(i);
                    partLUND_p[i] = sqrt(pow(partLUND_px[i], 2) + pow(partLUND_py[i], 2) + pow(partLUND_pz[i], 2));
                    partLUND_E[i] = std::sqrt(pow(partLUND_p[i], 2) + pow(partLUND_mass[i], 2));
                    partLUND_theta[i] = acos(partLUND_pz[i] / partLUND_p[i]);
                    partLUND_phi[i] = atan2(partLUND_py[i], partLUND_px[i]);
                }
            }
        }
        
        */
        
    }


    /// ////////////////////////////////////////////
}

/// //////////////////////////////////////////////////////////////////////////////////////////////////////
/// particle selection:

void select_electron(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:

                FD_eid_default_PID_check[i] = ele_default_PID_cut(i);
                FD_eid_charge_check[i] = ele_charge_cut(i);
                FD_eid_CC_nphe_check[i] = CC_nphe_cut(i);
                FD_eid_EC_outer_vs_EC_inner_check[i] = EC_outer_vs_EC_inner_cut(i);
                FD_eid_EC_sampling_fraction_check[i] = EC_sampling_fraction_cut(i);
                FD_eid_EC_hit_position_fiducial_check[i] = EC_hit_position_fiducial_cut_homogeneous(i);
                FD_eid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_eid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_eid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);
                FD_eid_DC_z_vertex_check[i] = DC_z_vertex_cut(i);

                // count statistics

                if (FD_eid_default_PID_check[i] == true)
                    FD_eid_default_PID_pass += 1;
                if (FD_eid_charge_check[i] == true)
                    FD_eid_charge_pass += 1;
                if (FD_eid_CC_nphe_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_CC_nphe_pass += 1;
                if (FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_EC_outer_vs_EC_inner_pass += 1;
                if (FD_eid_default_PID_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_EC_sampling_fraction_pass += 1;
                if (FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_EC_hit_position_fiducial_pass += 1;
                if (FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_eid_DC_z_vertex_check[i] && FD_eid_charge_check[i] && FD_eid_default_PID_check[i])
                    FD_eid_DC_z_vertex_pass += 1;

                if (FD_eid_default_PID_check[i] && FD_eid_charge_check[i] && FD_eid_EC_outer_vs_EC_inner_check[i] && FD_eid_EC_sampling_fraction_check[i] && FD_eid_EC_hit_position_fiducial_check[i] && FD_eid_DC_hit_position_region1_fiducial_check[i] && FD_eid_DC_hit_position_region2_fiducial_check[i] && FD_eid_DC_hit_position_region3_fiducial_check[i] && FD_eid_DC_z_vertex_check[i])
                {
                    FD_eid_all_pass += 1;
                    FD_eid_all_check[i] = true;
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward tagger:

            if (abs(part_status[i]) >= 1000 && abs(part_status[i]) < 2000)
            {
                FT_eid_charge_check[i] = FT_eid_charge_cut(i);
                FT_eid_PID_check[i] = FT_eid_PID_cut(i);
                FT_eid_FTCAL_fiducial_check[i] = FT_eid_FTCAL_fiducial_cut(i);
                FT_eid_FTTRK_fiducial_check[i] = FT_eid_FTTRK_fiducial_cut(i);
                FT_eid_FTHODO_fiducial_check[i] = FT_eid_FTHODO_fiducial_cut(i);
                FT_eid_energy_vs_radius_check[i] = FT_eid_energy_vs_radius_cut(i);

                if (FT_eid_charge_check[i])
                    FT_eid_charge_pass += 1;
                if (FT_eid_PID_check[i])
                    FT_eid_PID_pass += 1;
                if (FT_eid_FTCAL_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i])
                    FT_eid_FTCAL_fiducial_pass += 1;
                if (FT_eid_FTTRK_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i])
                    FT_eid_FTTRK_fiducial_pass += 1;
                if (FT_eid_FTHODO_fiducial_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i])
                    FT_eid_FTHODO_fiducial_pass += 1;
                if (FT_eid_energy_vs_radius_check[i] && FT_eid_charge_check[i] && FT_eid_PID_check[i])
                    FT_eid_energy_vs_radius_pass += 1;

                if (FT_eid_charge_check[i] && FT_eid_PID_check[i] && FT_eid_FTCAL_fiducial_check[i] && FT_eid_FTTRK_fiducial_check[i] && FT_eid_FTHODO_fiducial_check[i] && FT_eid_energy_vs_radius_check[i])
                {
                    FT_eid_all_pass += 1;
                    FT_eid_all_check[i] = true;
                }
            }

            // create vector with electron lorentz vectors if electron PID is succesfull:

            bool selector;

            if (use_FT == false)
            {
                FT_eid_all_check[i] = false;
            }
            if (use_FD == false)
            {
                FD_eid_all_check[i] = false;
            }

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            mom = 0;
            check = 0;

            for (int k = 0; k < Npart; k++)
            {
                if (k < BUFFER)
                {

                    if (use_own_PID_electron == true)
                    {
                        selector = (FD_eid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (FT_eid_all_check[k] && abs(part_status[i]) >= 1000 && abs(part_status[i]) < 2000);
                    }
                    else
                    {
                        selector = (FD_eid_default_PID_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (FT_eid_PID_check[k] && abs(part_status[i]) >= 1000 && abs(part_status[i]) < 2000);
                    }

                    if (selector)
                    {
                        for (int j = 0; j < i; j++)
                        {
                            if (k == e_ind[j])
                                check = -1;
                        }

                        if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                        {

                            mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                            e_ind[i] = k;
                            e_count += 1;
                        }
                        check = 0;
                    }
                }
            }
        }
    }

    /// ///////////////////////////////////////////////////////////////////
    /// assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (e_ind[i] != -1)
        {
            e_vx[i] = part_vx[e_ind[i]];
            e_vy[i] = part_vy[e_ind[i]];
            e_vz[i] = part_vz[e_ind[i]];
            e_beta[i] = part_beta[e_ind[i]];
            double p = sqrt(part_px[e_ind[i]] * part_px[e_ind[i]] + part_py[e_ind[i]] * part_py[e_ind[i]] + part_pz[e_ind[i]] * part_pz[e_ind[i]]);
            p4_ele[i].SetPxPyPzE(part_px[e_ind[i]], part_py[e_ind[i]], part_pz[e_ind[i]], sqrt(p * p + m_e * m_e));
            e_FTOF_sec[i] = part_FTOF_sector_layer2[e_ind[i]];
            e_PCAL_sec[i] = part_Cal_PCAL_sector[e_ind[i]];
            if (abs(part_status[e_ind[i]]) >= 1000 && abs(part_status[e_ind[i]]) < 2000)
                ele_detect[i] = 1;
            if (abs(part_status[e_ind[i]]) >= 2000 && abs(part_status[e_ind[i]]) < 4000)
                ele_detect[i] = 2;
            ele_chi2pid[i] = part_chi2pid[e_ind[i]];
            ele_dcx1[i] = part_DC_c1x[e_ind[i]];
            ele_dcy1[i] = part_DC_c1y[e_ind[i]];
            ele_dcx2[i] = part_DC_c2x[e_ind[i]];
            ele_dcy2[i] = part_DC_c2y[e_ind[i]];
            ele_dcx3[i] = part_DC_c3x[e_ind[i]];
            ele_dcy3[i] = part_DC_c3y[e_ind[i]];
            ele_dcedge1[i] = part_DC_edge1[e_ind[i]];
            ele_dcedge2[i] = part_DC_edge2[e_ind[i]];
            ele_dcedge3[i] = part_DC_edge3[e_ind[i]];
            ele_pcalv[i] = part_Cal_PCAL_lv[e_ind[i]];
            ele_pcalw[i] = part_Cal_PCAL_lw[e_ind[i]];         
            if (part_Cal_PCAL_sector[e_ind[i]] == 1)
                ele_sector[i] = 1;
            if (part_Cal_PCAL_sector[e_ind[i]] == 2)
                ele_sector[i] = 2;
            if (part_Cal_PCAL_sector[e_ind[i]] == 3)
                ele_sector[i] = 3;
            if (part_Cal_PCAL_sector[e_ind[i]] == 4)
                ele_sector[i] = 4;
            if (part_Cal_PCAL_sector[e_ind[i]] == 5)
                ele_sector[i] = 5;
            if (part_Cal_PCAL_sector[e_ind[i]] == 6)
                ele_sector[i] = 6;
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// proton selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_proton(int run)
{

    float mom = 0;
    float mom_new = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:

                FD_protid_default_PID_check[i] = prot_default_PID_cut(i);
                FD_protid_charge_check[i] = prot_charge_cut(i);  
                FD_protid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_protid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_protid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);
                FD_protid_delta_vz_check[i] = prot_delta_vz_cut(i);

                if (FD_protid_default_PID_check[i] == true)
                    FD_protid_default_PID_pass += 1;
                if (FD_protid_charge_check[i] == true)
                    FD_protid_charge_pass += 1;
                if (FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i])
                    FD_protid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i])
                    FD_protid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i])
                    FD_protid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_protid_delta_vz_check[i] && FD_protid_charge_check[i] && FD_protid_default_PID_check[i])
                    FD_protid_delta_vz_pass += 1;

                if (cut_maximum_probability_charged == true)
                {
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_maximum_probability_check[i] && FD_protid_delta_vz_check[i])
                    {
                        FD_protid_all_pass += 1;
                        FD_protid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true)
                {
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i] && FD_protid_delta_vz_check[i])
                    {
                        FD_protid_all_pass += 1;
                        FD_protid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true && cut_deltabeta_charged == true && cut_tofmass_charged == true)
                {
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_beta_check[i] && FD_protid_tofmass_check[i] && FD_protid_delta_beta_check[i] && FD_protid_delta_vz_check[i])
                    {
                        FD_protid_all_pass += 1;
                        FD_protid_all_check[i] = true;
                    }
                }
                else if (use_fiducial_charged == true)
                {
                    if (FD_protid_default_PID_check[i] && FD_protid_charge_check[i] && FD_protid_DC_hit_position_region1_fiducial_check[i] && FD_protid_DC_hit_position_region2_fiducial_check[i] && FD_protid_DC_hit_position_region3_fiducial_check[i] && FD_protid_delta_vz_check[i])
                    {
                        FD_protid_all_pass += 1;
                        FD_protid_all_check[i] = true;
                    }
                }
                else
                {
                    if (FD_protid_default_PID_check[i])
                    {
                        FD_protid_all_pass += 1;
                        FD_protid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                CD_protid_default_PID_check[i] = prot_default_PID_cut(i);
                CD_protid_charge_check[i] = prot_charge_cut(i);
                CD_protid_delta_vz_check[i] = CD_prot_delta_vz_cut(i);

                if (CD_protid_default_PID_check[i])
                    CD_protid_default_PID_pass += 1;
                if (CD_protid_charge_check[i])
                    CD_protid_charge_pass += 1;
                if (CD_protid_delta_vz_check[i] && CD_protid_charge_check[i] && CD_protid_default_PID_check[i])
                    CD_protid_delta_vz_pass += 1;
                else if (use_CD_fiducial_charged == true)
                {
                    if (CD_protid_default_PID_check[i] && CD_protid_charge_check[i] && CD_protid_delta_vz_check[i])
                    {
                        CD_protid_all_pass += 1;
                        CD_protid_all_check[i] = true;
                    }
                }
                else
                {
                    if (CD_protid_default_PID_check[i])
                    {
                        CD_protid_all_pass += 1;
                        CD_protid_all_check[i] = true;
                    }
                }

                /// ///////////////////////////////////////////////////////////////////////
                /// 2 GeV condition:

                if (Ebeam < 3)
                {
                    if (CD_protid_charge_check[i])
                    {
                        CD_protid_all_check[i] = true;
                    }
                }

                /// //////////////////////////////////////////////////////////////////////
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create proton selector

            bool selector;

            if (use_FD == false)
                FD_protid_all_check[i] = false;
            if (use_CD == false)
                CD_protid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_protid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_protid_all_check[k] && part_status[i] >= 4000);

                        ///

                        if (selector)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == p_ind[j])
                                    check = -1;
                            }
                            mom_new = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                            if (mom_new > mom && mom_new > 0 && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                p_ind[i] = k;
                                p_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (p_ind[i] != -1)
        {
            p_vx[i] = part_vx[p_ind[i]];
            p_vy[i] = part_vy[p_ind[i]];
            p_vz[i] = part_vz[p_ind[i]];
            p_beta[i] = part_beta[p_ind[i]];
            double p = sqrt(part_px[p_ind[i]] * part_px[p_ind[i]] + part_py[p_ind[i]] * part_py[p_ind[i]] + part_pz[p_ind[i]] * part_pz[p_ind[i]]);
            p4_prot[i].SetPxPyPzE(part_px[p_ind[i]], part_py[p_ind[i]], part_pz[p_ind[i]], sqrt(p * p + m_p * m_p));
            p_FTOF_sec[i] = part_FTOF_sector_layer2[p_ind[i]];
            p_PCAL_sec[i] = part_Cal_PCAL_sector[p_ind[i]];
            if (abs(part_status[p_ind[i]]) >= 2000 && abs(part_status[p_ind[i]]) < 4000)
                prot_detect[i] = 2;
            if (part_status[p_ind[i]] >= 4000)
                prot_detect[i] = 3;
            prot_chi2pid[i] = part_chi2pid[p_ind[i]];
            prot_dcx1[i] = part_DC_c1x[p_ind[i]];
            prot_dcy1[i] = part_DC_c1y[p_ind[i]];
            prot_dcz1[i] = part_DC_c1z[p_ind[i]];
            prot_dcedge1[i] = part_DC_edge1[p_ind[i]];
            prot_dcedge2[i] = part_DC_edge2[p_ind[i]];
            prot_dcedge3[i] = part_DC_edge3[p_ind[i]];
            if (part_DC_sector[p_ind[i]] == 1)
                prot_sector[i] = 1;
            if (part_DC_sector[p_ind[i]] == 2)
                prot_sector[i] = 2;
            if (part_DC_sector[p_ind[i]] == 3)
                prot_sector[i] = 3;
            if (part_DC_sector[p_ind[i]] == 4)
                prot_sector[i] = 4;
            if (part_DC_sector[p_ind[i]] == 5)
                prot_sector[i] = 5;
            if (part_DC_sector[p_ind[i]] == 6)
                prot_sector[i] = 6;
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// neutron selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_neutron(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:

                FD_neutrid_default_PID_check[i] = neutr_default_PID_cut(i);
                FD_neutrid_charge_check[i] = neutr_charge_cut(i);
                FD_neutrid_beta_check[i] = neutr_beta_cut(i, run);
                // FD_neutrid_delta_beta_check[i] = neutr_delta_beta_cut(i, run);
                // FD_neutrid_tofmass_check[i] = neutr_tofmass_cut(i, run);
                FD_neutrid_delta_vz_check[i] = neutr_delta_vz_cut(i);

                if (FD_neutrid_default_PID_check[i] == true)
                    FD_neutrid_default_PID_pass += 1;
                if (FD_neutrid_charge_check[i] == true)
                    FD_neutrid_charge_pass += 1;
                if (FD_neutrid_beta_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i])
                    FD_neutrid_beta_pass += 1;
                // if(FD_neutrid_delta_beta_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 		FD_neutrid_delta_beta_pass += 1;
                // if(FD_neutrid_tofmass_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i]) 		FD_neutrid_tofmass_pass += 1;
                if (FD_neutrid_delta_vz_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_default_PID_check[i])
                    FD_neutrid_delta_vz_pass += 1;

                if (cut_extend_neutron == true)
                {
                    if (part_p[i] > 0.05 && part_p[i] < 2.0 && FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i] && FD_neutrid_delta_vz_check[i])
                    {
                        FD_neutrid_all_pass += 1;
                        FD_neutrid_all_check[i] = true;
                    }
                    if (part_p[i] >= 2.0 && FD_neutrid_charge_check[i] && FD_neutrid_delta_vz_check[i])
                    {
                        FD_neutrid_all_pass += 1;
                        FD_neutrid_all_check[i] = true;
                    }
                }
                else
                {
                    if (part_p[i] > 0.05 && FD_neutrid_default_PID_check[i] && FD_neutrid_charge_check[i] && FD_neutrid_beta_check[i] && FD_neutrid_delta_vz_check[i])
                    {
                        FD_neutrid_all_pass += 1;
                        FD_neutrid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                CD_neutrid_default_PID_check[i] = neutr_default_PID_cut(i);
                CD_neutrid_charge_check[i] = neutr_charge_cut(i);
                CD_neutrid_beta_check[i] = CD_neutr_beta_cut(i, run);
                CD_neutrid_delta_vz_check[i] = CD_neutr_delta_vz_cut(i);

                if (CD_neutrid_default_PID_check[i])
                    CD_neutrid_default_PID_pass += 1;
                if (CD_neutrid_charge_check[i])
                    CD_neutrid_charge_pass += 1;
                if (CD_neutrid_beta_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_default_PID_check[i])
                    CD_neutrid_beta_pass += 1;
                if (CD_neutrid_delta_vz_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_default_PID_check[i])
                    CD_neutrid_delta_vz_pass += 1;

                if (part_p[i] > 0.03 && CD_neutrid_default_PID_check[i] && CD_neutrid_charge_check[i] && CD_neutrid_beta_check[i] && CD_neutrid_delta_vz_check[i])
                {
                    CD_neutrid_all_pass += 1;
                    CD_neutrid_all_check[i] = true;
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create neutron selector

            bool selector;

            if (use_FD == false)
                FD_neutrid_all_check[i] = false;
            if (use_CD == false)
                CD_neutrid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_neutrid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_neutrid_all_check[k] && part_status[i] >= 4000);

                        if (selector)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == n_ind[j])
                                    check = -1;
                            }
                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                n_ind[i] = k;
                                n_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (n_ind[i] != -1)
        {
            n_vx[i] = part_vx[n_ind[i]];
            n_vy[i] = part_vy[n_ind[i]];
            n_vz[i] = part_vz[n_ind[i]];
            n_beta[i] = part_beta[n_ind[i]];
            double p = sqrt(part_px[n_ind[i]] * part_px[n_ind[i]] + part_py[n_ind[i]] * part_py[n_ind[i]] + part_pz[n_ind[i]] * part_pz[n_ind[i]]);
            p4_neutr[i].SetPxPyPzE(part_px[n_ind[i]], part_py[n_ind[i]], part_pz[n_ind[i]], sqrt(p * p + m_n * m_n));
            double cx = part_px[n_ind[i]] / p;
            double cy = part_py[n_ind[i]] / p;
            double cz = part_pz[n_ind[i]] / p;
            neutr_chi2pid[i] = part_chi2pid[n_ind[i]];
            if (p > 2.0)
            {
                double beta_neutr = part_beta[n_ind[i]];
                double gamma_neutr = sqrt(1 / (1 - beta_neutr * beta_neutr));
                double E = gamma_neutr * m_n;
                p = sqrt(E * E - m_n * m_n);
                p4_neutr[i].SetPxPyPzE(p * cx, p * cy, p * cz, E);
            }
            if (abs(part_status[n_ind[i]]) >= 2000 && abs(part_status[n_ind[i]]) < 4000)
                neutr_detect[i] = 2;
            if (part_status[n_ind[i]] >= 4000)
                neutr_detect[i] = 3;
            if (part_Cal_PCAL_sector[n_ind[i]] == 1)
                neutr_sector[i] = 1;
            if (part_Cal_PCAL_sector[n_ind[i]] == 2)
                neutr_sector[i] = 2;
            if (part_Cal_PCAL_sector[n_ind[i]] == 3)
                neutr_sector[i] = 3;
            if (part_Cal_PCAL_sector[n_ind[i]] == 4)
                neutr_sector[i] = 4;
            if (part_Cal_PCAL_sector[n_ind[i]] == 5)
                neutr_sector[i] = 5;
            if (part_Cal_PCAL_sector[n_ind[i]] == 6)
                neutr_sector[i] = 6;
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// pip selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_pip(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:
                FD_pipid_default_PID_check[i] = pip_default_PID_cut(i);
                FD_pipid_charge_check[i] = pip_charge_cut(i); 
                FD_pipid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_pipid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_pipid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);
                FD_pipid_delta_vz_check[i] = pip_delta_vz_cut(i);

                if (FD_pipid_default_PID_check[i] == true)
                    FD_pipid_default_PID_pass += 1;
                if (FD_pipid_charge_check[i] == true)
                    FD_pipid_charge_pass += 1;
                if (FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_pipid_beta_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_beta_pass += 1;
                if (FD_pipid_delta_beta_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_delta_beta_pass += 1;
                if (FD_pipid_tofmass_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_tofmass_pass += 1;
                if (FD_pipid_maximum_probability_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_maximum_probability_pass += 1;
                if (FD_pipid_delta_vz_check[i] && FD_pipid_charge_check[i] && FD_pipid_default_PID_check[i])
                    FD_pipid_delta_vz_pass += 1;

                if (cut_maximum_probability_charged == true)
                {
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_maximum_probability_check[i] && FD_pipid_delta_vz_check[i])
                    {
                        FD_pipid_all_pass += 1;
                        FD_pipid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true)
                {
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i] && FD_pipid_delta_vz_check[i])
                    {
                        FD_pipid_all_pass += 1;
                        FD_pipid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true && cut_deltabeta_charged == true && cut_tofmass_charged == true)
                {
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_beta_check[i] && FD_pipid_delta_beta_check[i] && FD_pipid_tofmass_check[i] && FD_pipid_delta_vz_check[i])
                    {
                        FD_pipid_all_pass += 1;
                        FD_pipid_all_check[i] = true;
                    }
                }
                else if (use_fiducial_charged == true)
                {
                    if (FD_pipid_default_PID_check[i] && FD_pipid_charge_check[i] && FD_pipid_DC_hit_position_region1_fiducial_check[i] && FD_pipid_DC_hit_position_region2_fiducial_check[i] && FD_pipid_DC_hit_position_region3_fiducial_check[i] && FD_pipid_delta_vz_check[i])
                    {
                        FD_pipid_all_pass += 1;
                        FD_pipid_all_check[i] = true;
                    }
                }
                else
                {
                    if (FD_pipid_default_PID_check[i])
                    {
                        FD_pipid_all_pass += 1;
                        FD_pipid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                CD_pipid_default_PID_check[i] = pip_default_PID_cut(i);
                CD_pipid_charge_check[i] = pip_charge_cut(i);
                CD_pipid_delta_vz_check[i] = CD_pip_delta_vz_cut(i);

                if (CD_pipid_default_PID_check[i])
                    CD_pipid_default_PID_pass += 1;
                if (CD_pipid_charge_check[i])
                    CD_pipid_charge_pass += 1;
                if (CD_pipid_beta_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i])
                    CD_pipid_beta_pass += 1;
                if (CD_pipid_maximum_probability_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i])
                    CD_pipid_maximum_probability_pass += 1;
                if (CD_pipid_delta_vz_check[i] && CD_pipid_charge_check[i] && CD_pipid_default_PID_check[i])
                    CD_pipid_delta_vz_pass += 1;

                if (CD_cut_maximum_probability_charged == true)
                {
                    if (CD_pipid_default_PID_check[i] && CD_pipid_charge_check[i] && CD_pipid_maximum_probability_check[i] && CD_pipid_delta_vz_check[i])
                    {
                        CD_pipid_all_pass += 1;
                        CD_pipid_all_check[i] = true;
                    }
                }
                else if (CD_cut_beta_vs_p_charged == true)
                {
                    if (CD_pipid_default_PID_check[i] && CD_pipid_charge_check[i] && CD_pipid_beta_check[i] && CD_pipid_delta_vz_check[i])
                    {
                        CD_pipid_all_pass += 1;
                        CD_pipid_all_check[i] = true;
                    }
                }
                else if (use_CD_fiducial_charged == true)
                {
                    if (CD_pipid_default_PID_check[i] && CD_pipid_charge_check[i] && CD_pipid_delta_vz_check[i])
                    {
                        CD_pipid_all_pass += 1;
                        CD_pipid_all_check[i] = true;
                    }
                }
                else
                {
                    if (CD_pipid_default_PID_check[i])
                    {
                        CD_pipid_all_pass += 1;
                        CD_pipid_all_check[i] = true;
                    }
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create pip selector

            bool selector;

            if (use_FD == false)
                FD_pipid_all_check[i] = false;
            if (use_CD == false)
                CD_pipid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_pipid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_pipid_all_check[k] && part_status[i] >= 4000);

                        if (selector == true)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == pip_ind[j])
                                    check = -1;
                            }
                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                pip_ind[i] = k;
                                pip_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (pip_ind[i] != -1)
        {
            pip_vx[i] = part_vx[pip_ind[i]];
            pip_vy[i] = part_vy[pip_ind[i]];
            pip_vz[i] = part_vz[pip_ind[i]];
            pip_beta[i] = part_beta[pip_ind[i]];
            double p = sqrt(part_px[pip_ind[i]] * part_px[pip_ind[i]] + part_py[pip_ind[i]] * part_py[pip_ind[i]] + part_pz[pip_ind[i]] * part_pz[pip_ind[i]]);
            p4_pip[i].SetPxPyPzE(part_px[pip_ind[i]], part_py[pip_ind[i]], part_pz[pip_ind[i]], sqrt(p * p + m_pip * m_pip));
            pip_FTOF_sec[i] = part_FTOF_sector_layer2[pip_ind[i]];
            if (abs(part_status[pip_ind[i]]) >= 2000 && abs(part_status[pip_ind[i]]) < 4000)
                pip_detect[i] = 2;
            if (part_status[pip_ind[i]] >= 4000)
                pip_detect[i] = 3;
            pip_chi2pid[i] = part_chi2pid[pip_ind[i]];
            pip_dcx1[i] = part_DC_c1x[pip_ind[i]];
            pip_dcy1[i] = part_DC_c1y[pip_ind[i]];
            pip_dcx2[i] = part_DC_c2x[pip_ind[i]];
            pip_dcy2[i] = part_DC_c2y[pip_ind[i]];
            pip_dcx3[i] = part_DC_c3x[pip_ind[i]];
            pip_dcy3[i] = part_DC_c3y[pip_ind[i]];
            pip_dcedge1[i] = part_DC_edge1[pip_ind[i]];
            pip_dcedge2[i] = part_DC_edge2[pip_ind[i]];
            pip_dcedge3[i] = part_DC_edge3[pip_ind[i]];
            if(part_DC_sector[pip_ind[i]] == 1) pip_sector[i] = 1;
            if(part_DC_sector[pip_ind[i]] == 2) pip_sector[i] = 2; 
            if(part_DC_sector[pip_ind[i]] == 3) pip_sector[i] = 3; 
            if(part_DC_sector[pip_ind[i]] == 4) pip_sector[i] = 4; 
            if(part_DC_sector[pip_ind[i]] == 5) pip_sector[i] = 5; 
            if(part_DC_sector[pip_ind[i]] == 6) pip_sector[i] = 6;
            // pip_RICHbest[i] = part_RICH_best[pip_ind[i]];
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// pim selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_pim(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:

                FD_pimid_default_PID_check[i] = pim_default_PID_cut(i);
                FD_pimid_charge_check[i] = pim_charge_cut(i);
                FD_pimid_ele_reject_check[i] = pim_ele_reject_cut(i);
                FD_pimid_EC_outer_vs_EC_inner_check[i] = pim_EC_outer_vs_EC_inner_cut(i);
                FD_pimid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_pimid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_pimid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);  
                FD_pimid_delta_vz_check[i] = pim_delta_vz_cut(i);

                if (FD_pimid_default_PID_check[i] == true)
                    FD_pimid_default_PID_pass += 1;
                if (FD_pimid_charge_check[i] == true)
                    FD_pimid_charge_pass += 1;
                if (FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_pimid_beta_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_beta_pass += 1;
                if (FD_pimid_delta_beta_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_delta_beta_pass += 1;
                if (FD_pimid_tofmass_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_tofmass_pass += 1;
                if (FD_pimid_maximum_probability_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_maximum_probability_pass += 1;
                if (FD_pimid_delta_vz_check[i] && FD_pimid_charge_check[i] && FD_pimid_default_PID_check[i])
                    FD_pimid_delta_vz_pass += 1;

                if (cut_maximum_probability_charged == true)
                {
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_maximum_probability_check[i] && FD_pimid_delta_vz_check[i])
                    {
                        FD_pimid_all_pass += 1;
                        FD_pimid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true)
                {
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_beta_check[i] && FD_pimid_delta_vz_check[i])
                    {
                        FD_pimid_all_pass += 1;
                        FD_pimid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true && cut_deltabeta_charged == true && cut_tofmass_charged == true)
                {
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_beta_check[i] && FD_pimid_delta_beta_check[i] && FD_pimid_tofmass_check[i] && FD_pimid_delta_vz_check[i])
                    {
                        FD_pimid_all_pass += 1;
                        FD_pimid_all_check[i] = true;
                    }
                }
                else if (use_fiducial_charged == true)
                {
                    if (FD_pimid_default_PID_check[i] && FD_pimid_charge_check[i] && FD_pimid_DC_hit_position_region1_fiducial_check[i] && FD_pimid_DC_hit_position_region2_fiducial_check[i] && FD_pimid_DC_hit_position_region3_fiducial_check[i] && FD_pimid_delta_vz_check[i])
                    {
                        FD_pimid_all_pass += 1;
                        FD_pimid_all_check[i] = true;
                    }
                }
                else
                {
                    if (FD_pimid_default_PID_check[i])
                    {
                        FD_pimid_all_pass += 1;
                        FD_pimid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                double beta_charge_central = Beta_charged_central(i, run);

                CD_pimid_default_PID_check[i] = pim_default_PID_cut(i);
                CD_pimid_charge_check[i] = pim_charge_cut(i);
                CD_pimid_delta_vz_check[i] = CD_pim_delta_vz_cut(i);

                if (CD_pimid_default_PID_check[i])
                    CD_pimid_default_PID_pass += 1;
                if (CD_pimid_charge_check[i])
                    CD_pimid_charge_pass += 1;
                if (CD_pimid_beta_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i])
                    CD_pimid_beta_pass += 1;
                if (CD_pimid_maximum_probability_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i])
                    CD_pimid_maximum_probability_pass += 1;
                if (CD_pimid_delta_vz_check[i] && CD_pimid_charge_check[i] && CD_pimid_default_PID_check[i])
                    CD_pimid_delta_vz_pass += 1;

                if (CD_cut_maximum_probability_charged == true)
                {
                    if (CD_pimid_default_PID_check[i] && CD_pimid_charge_check[i] && CD_pimid_maximum_probability_check[i] && CD_pimid_delta_vz_check[i])
                    {
                        CD_pimid_all_pass += 1;
                        CD_pimid_all_check[i] = true;
                    }
                }
                else if (CD_cut_beta_vs_p_charged == true)
                {
                    if (CD_pimid_default_PID_check[i] && CD_pimid_charge_check[i] && CD_pimid_beta_check[i] && CD_pimid_delta_vz_check[i])
                    {
                        CD_pimid_all_pass += 1;
                        CD_pimid_all_check[i] = true;
                    }
                }
                else if (use_CD_fiducial_charged == true)
                {
                    if (CD_pimid_default_PID_check[i] && CD_pimid_charge_check[i] && CD_pimid_delta_vz_check[i])
                    {
                        CD_pimid_all_pass += 1;
                        CD_pimid_all_check[i] = true;
                    }
                }
                else
                {
                    if (CD_pimid_default_PID_check[i])
                    {
                        CD_pimid_all_pass += 1;
                        CD_pimid_all_check[i] = true;
                    }
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create pim selector

            bool selector;

            if (use_FD == false)
                FD_pimid_all_check[i] = false;
            if (use_CD == false)
                CD_pimid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_pimid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_pimid_all_check[k] && part_status[i] >= 4000);

                        if (selector)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == pim_ind[j])
                                    check = -1;
                            }
                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                pim_ind[i] = k;
                                pim_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (pim_ind[i] != -1)
        {
            pim_vx[i] = part_vx[pim_ind[i]];
            pim_vy[i] = part_vy[pim_ind[i]];
            pim_vz[i] = part_vz[pim_ind[i]];
            pim_beta[i] = part_beta[pim_ind[i]];
            double p = sqrt(part_px[pim_ind[i]] * part_px[pim_ind[i]] + part_py[pim_ind[i]] * part_py[pim_ind[i]] + part_pz[pim_ind[i]] * part_pz[pim_ind[i]]);
            p4_pim[i].SetPxPyPzE(part_px[pim_ind[i]], part_py[pim_ind[i]], part_pz[pim_ind[i]], sqrt(p * p + m_pim * m_pim));
            pim_FTOF_sec[i] = part_FTOF_sector_layer2[pim_ind[i]];
            if (abs(part_status[pim_ind[i]]) >= 2000 && abs(part_status[pim_ind[i]]) < 4000)
                pim_detect[i] = 2;
            if (part_status[pim_ind[i]] >= 4000)
                pim_detect[i] = 3;
            pim_chi2pid[i] = part_chi2pid[pim_ind[i]];
            pim_dcx1[i] = part_DC_c1x[pim_ind[i]];
            pim_dcy1[i] = part_DC_c1y[pim_ind[i]];
            pim_dcx2[i] = part_DC_c2x[pim_ind[i]];
            pim_dcy2[i] = part_DC_c2y[pim_ind[i]];
            pim_dcx3[i] = part_DC_c3x[pim_ind[i]];
            pim_dcy3[i] = part_DC_c3y[pim_ind[i]];
            pim_dcedge1[i] = part_DC_edge1[pim_ind[i]];
            pim_dcedge2[i] = part_DC_edge2[pim_ind[i]];
            pim_dcedge3[i] = part_DC_edge3[pim_ind[i]];
            if (part_DC_sector[pim_ind[i]] == 1)
                pim_sector[i] = 1;
            if (part_DC_sector[pim_ind[i]] == 2)
                pim_sector[i] = 2;
            if (part_DC_sector[pim_ind[i]] == 3)
                pim_sector[i] = 3;
            if (part_DC_sector[pim_ind[i]] == 4)
                pim_sector[i] = 4;
            if (part_DC_sector[pim_ind[i]] == 5)
                pim_sector[i] = 5;
            if (part_DC_sector[pim_ind[i]] == 6)
                pim_sector[i] = 6;
            // pim_RICHbest[i] = part_RICH_best[pim_ind[i]];
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Kp selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_Kplus(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000 && part_DC_sector[i] != 0)
            {

                // PID checks:

                FD_Kpid_default_PID_check[i] = Kp_default_PID_cut(i);
                FD_Kpid_charge_check[i] = Kp_charge_cut(i);
                FD_Kpid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_Kpid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_Kpid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);
                FD_Kpid_delta_vz_check[i] = Kp_delta_vz_cut(i);

                if (FD_Kpid_default_PID_check[i] == true)
                    FD_Kpid_default_PID_pass += 1;
                if (FD_Kpid_charge_check[i] == true)
                    FD_Kpid_charge_pass += 1;
                if (FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_Kpid_beta_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_beta_pass += 1;
                if (FD_Kpid_delta_beta_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_delta_beta_pass += 1;
                if (FD_Kpid_tofmass_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_tofmass_pass += 1;
                if (FD_Kpid_maximum_probability_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_maximum_probability_pass += 1;
                if (FD_Kpid_delta_vz_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_default_PID_check[i])
                    FD_Kpid_delta_vz_pass += 1;

                if (cut_maximum_probability_charged == true)
                {
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_maximum_probability_check[i] && FD_Kpid_delta_vz_check[i])
                    {
                        FD_Kpid_all_pass += 1;
                        FD_Kpid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true)
                {
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i] && FD_Kpid_delta_vz_check[i])
                    {
                        FD_Kpid_all_pass += 1;
                        FD_Kpid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true && cut_deltabeta_charged == true && cut_tofmass_charged == true)
                {
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_beta_check[i] && FD_Kpid_delta_beta_check[i] && FD_Kpid_tofmass_check[i] && FD_Kpid_delta_vz_check[i])
                    {
                        FD_Kpid_all_pass += 1;
                        FD_Kpid_all_check[i] = true;
                    }
                }
                else if (use_fiducial_charged == true)
                {
                    if (FD_Kpid_default_PID_check[i] && FD_Kpid_charge_check[i] && FD_Kpid_DC_hit_position_region1_fiducial_check[i] && FD_Kpid_DC_hit_position_region2_fiducial_check[i] && FD_Kpid_DC_hit_position_region3_fiducial_check[i] && FD_Kpid_delta_vz_check[i])
                    {
                        FD_Kpid_all_pass += 1;
                        FD_Kpid_all_check[i] = true;
                    }
                }
                else
                {
                    if (FD_Kpid_default_PID_check[i])
                    {
                        FD_Kpid_all_pass += 1;
                        FD_Kpid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                CD_Kpid_default_PID_check[i] = Kp_default_PID_cut(i);
                CD_Kpid_charge_check[i] = Kp_charge_cut(i);
                CD_Kpid_delta_vz_check[i] = CD_Kp_delta_vz_cut(i);

                if (CD_Kpid_default_PID_check[i])
                    CD_Kpid_default_PID_pass += 1;
                if (CD_Kpid_charge_check[i])
                    CD_Kpid_charge_pass += 1;
                if (CD_Kpid_beta_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i])
                    CD_Kpid_beta_pass += 1;
                if (CD_Kpid_maximum_probability_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i])
                    CD_Kpid_maximum_probability_pass += 1;
                if (CD_Kpid_delta_vz_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_default_PID_check[i])
                    CD_Kpid_delta_vz_pass += 1;

                if (CD_cut_maximum_probability_charged == true)
                {
                    if (CD_Kpid_default_PID_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_maximum_probability_check[i] && CD_Kpid_delta_vz_check[i])
                    {
                        CD_Kpid_all_pass += 1;
                        CD_Kpid_all_check[i] = true;
                    }
                }
                else if (CD_cut_beta_vs_p_charged == true)
                {
                    if (CD_Kpid_default_PID_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_beta_check[i] && CD_Kpid_delta_vz_check[i])
                    {
                        CD_Kpid_all_pass += 1;
                        CD_Kpid_all_check[i] = true;
                    }
                }
                else if (use_CD_fiducial_charged == true)
                {
                    if (CD_Kpid_default_PID_check[i] && CD_Kpid_charge_check[i] && CD_Kpid_delta_vz_check[i])
                    {
                        CD_Kpid_all_pass += 1;
                        CD_Kpid_all_check[i] = true;
                    }
                }
                else
                {
                    if (CD_Kpid_default_PID_check[i])
                    {
                        CD_Kpid_all_pass += 1;
                        CD_Kpid_all_check[i] = true;
                    }
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create Kp selector

            bool selector;

            if (use_FD == false)
                FD_Kpid_all_check[i] = false;
            if (use_CD == false)
                CD_Kpid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_Kpid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_Kpid_all_check[k] && part_status[i] >= 4000);

                        if (selector)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == Kp_ind[j])
                                    check = -1;
                            }
                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                Kp_ind[i] = k;
                                Kp_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (Kp_ind[i] != -1)
        {
            Kp_vx[i] = part_vx[Kp_ind[i]];
            Kp_vy[i] = part_vy[Kp_ind[i]];
            Kp_vz[i] = part_vz[Kp_ind[i]];
            Kp_beta[i] = part_beta[Kp_ind[i]];
            double p = sqrt(part_px[Kp_ind[i]] * part_px[Kp_ind[i]] + part_py[Kp_ind[i]] * part_py[Kp_ind[i]] + part_pz[Kp_ind[i]] * part_pz[Kp_ind[i]]);
            p4_Kp[i].SetPxPyPzE(part_px[Kp_ind[i]], part_py[Kp_ind[i]], part_pz[Kp_ind[i]], sqrt(p * p + m_Kp * m_Kp));
            Kp_FTOF_sec[i] = part_FTOF_sector_layer2[Kp_ind[i]];
            if (abs(part_status[Kp_ind[i]]) >= 2000 && abs(part_status[Kp_ind[i]]) < 4000)
                Kp_detect[i] = 2;
            if (part_status[Kp_ind[i]] >= 4000)
                Kp_detect[i] = 3;
            Kp_chi2pid[i] = part_chi2pid[Kp_ind[i]];
            if (part_DC_sector[Kp_ind[i]] == 1)
                Kp_sector[i] = 1;
            if (part_DC_sector[Kp_ind[i]] == 2)
                Kp_sector[i] = 2;
            if (part_DC_sector[Kp_ind[i]] == 3)
                Kp_sector[i] = 3;
            if (part_DC_sector[Kp_ind[i]] == 4)
                Kp_sector[i] = 4;
            if (part_DC_sector[Kp_ind[i]] == 5)
                Kp_sector[i] = 5;
            if (part_DC_sector[Kp_ind[i]] == 6)
                Kp_sector[i] = 6;
            // Kp_RICHbest[i] = part_RICH_best[Kp_ind[i]];
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Km selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_Kminus(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Forward detector:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                // PID checks:

                FD_Kmid_default_PID_check[i] = Km_default_PID_cut(i);
                FD_Kmid_charge_check[i] = Km_charge_cut(i);
                FD_Kmid_ele_reject_check[i] = Km_ele_reject_cut(i);
                FD_Kmid_EC_outer_vs_EC_inner_check[i] = Km_EC_outer_vs_EC_inner_cut(i);
                FD_Kmid_DC_hit_position_region1_fiducial_check[i] = DC_fiducial_cut_edge(i, 1);
                FD_Kmid_DC_hit_position_region2_fiducial_check[i] = DC_fiducial_cut_edge(i, 2);
                FD_Kmid_DC_hit_position_region3_fiducial_check[i] = DC_fiducial_cut_edge(i, 3);
                FD_Kmid_delta_vz_check[i] = Km_delta_vz_cut(i);

                if (FD_Kmid_default_PID_check[i] == true)
                    FD_Kmid_default_PID_pass += 1;
                if (FD_Kmid_charge_check[i] == true)
                    FD_Kmid_charge_pass += 1;
                if (FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_DC_hit_position_region1_fiducial_pass += 1;
                if (FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_DC_hit_position_region2_fiducial_pass += 1;
                if (FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_DC_hit_position_region3_fiducial_pass += 1;
                if (FD_Kmid_beta_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_beta_pass += 1;
                if (FD_Kmid_delta_beta_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_delta_beta_pass += 1;
                if (FD_Kmid_tofmass_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_tofmass_pass += 1;
                if (FD_Kmid_maximum_probability_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_maximum_probability_pass += 1;
                if (FD_Kmid_delta_vz_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_default_PID_check[i])
                    FD_Kmid_delta_vz_pass += 1;

                if (cut_maximum_probability_charged == true)
                {
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_maximum_probability_check[i] && FD_Kmid_delta_vz_check[i])
                    {
                        FD_Kmid_all_pass += 1;
                        FD_Kmid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true)
                {
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_beta_check[i] && FD_Kmid_delta_vz_check[i])
                    {
                        FD_Kmid_all_pass += 1;
                        FD_Kmid_all_check[i] = true;
                    }
                }
                else if (cut_beta_vs_p_charged == true && cut_deltabeta_charged == true && cut_tofmass_charged == true)
                {
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_beta_check[i] && FD_Kmid_delta_beta_check[i] && FD_Kmid_tofmass_check[i] && FD_Kmid_delta_vz_check[i])
                    {
                        FD_Kmid_all_pass += 1;
                        FD_Kmid_all_check[i] = true;
                    }
                }
                else if (use_fiducial_charged == true)
                {
                    if (FD_Kmid_default_PID_check[i] && FD_Kmid_charge_check[i] && FD_Kmid_DC_hit_position_region1_fiducial_check[i] && FD_Kmid_DC_hit_position_region2_fiducial_check[i] && FD_Kmid_DC_hit_position_region3_fiducial_check[i] && FD_Kmid_delta_vz_check[i])
                    {
                        FD_Kmid_all_pass += 1;
                        FD_Kmid_all_check[i] = true;
                    }
                }
                else
                {
                    if (FD_Kmid_default_PID_check[i])
                    {
                        FD_Kmid_all_pass += 1;
                        FD_Kmid_all_check[i] = true;
                    }
                }
            }

            /// ////////////////////////////////////////////////////////////////////////////////////////////
            /// Central detector:

            if (part_status[i] >= 4000)
            {

                // PID checks:

                double beta_charge_central = Beta_charged_central(i, run);

                CD_Kmid_default_PID_check[i] = Km_default_PID_cut(i);
                CD_Kmid_charge_check[i] = Km_charge_cut(i);
                CD_Kmid_delta_vz_check[i] = CD_Km_delta_vz_cut(i);

                if (CD_Kmid_default_PID_check[i])
                    CD_Kmid_default_PID_pass += 1;
                if (CD_Kmid_charge_check[i])
                    CD_Kmid_charge_pass += 1;
                if (CD_Kmid_beta_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i])
                    CD_Kmid_beta_pass += 1;
                if (CD_Kmid_maximum_probability_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i])
                    CD_Kmid_maximum_probability_pass += 1;
                if (CD_Kmid_delta_vz_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_default_PID_check[i])
                    CD_Kmid_delta_vz_pass += 1;

                if (CD_cut_maximum_probability_charged == true)
                {
                    if (CD_Kmid_default_PID_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_maximum_probability_check[i] && CD_Kmid_delta_vz_check[i])
                    {
                        CD_Kmid_all_pass += 1;
                        CD_Kmid_all_check[i] = true;
                    }
                }
                else if (CD_cut_beta_vs_p_charged == true)
                {
                    if (CD_Kmid_default_PID_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_beta_check[i] && CD_Kmid_delta_vz_check[i])
                    {
                        CD_Kmid_all_pass += 1;
                        CD_Kmid_all_check[i] = true;
                    }
                }
                else if (use_CD_fiducial_charged == true)
                {
                    if (CD_Kmid_default_PID_check[i] && CD_Kmid_charge_check[i] && CD_Kmid_delta_vz_check[i])
                    {
                        CD_Kmid_all_pass += 1;
                        CD_Kmid_all_check[i] = true;
                    }
                }
                else
                {
                    if (CD_Kmid_default_PID_check[i])
                    {
                        CD_Kmid_all_pass += 1;
                        CD_Kmid_all_check[i] = true;
                    }
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create Km selector

            bool selector;

            if (use_FD == false)
                FD_Kmid_all_check[i] = false;
            if (use_CD == false)
                CD_Kmid_all_check[i] = false;

            /// ////////////////////////////////////////////////////////////////
            /// pick particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_Kmid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (CD_Kmid_all_check[k] && part_status[i] >= 4000);

                        if (selector)
                        {
                            for (int j = 0; j < i; j++)
                            {
                                if (k == Km_ind[j])
                                    check = -1;
                            }
                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                Km_ind[i] = k;
                                Km_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (Km_ind[i] != -1)
        {
            Km_vx[i] = part_vx[Km_ind[i]];
            Km_vy[i] = part_vy[Km_ind[i]];
            Km_vz[i] = part_vz[Km_ind[i]];
            Km_beta[i] = part_beta[Km_ind[i]];
            double p = sqrt(part_px[Km_ind[i]] * part_px[Km_ind[i]] + part_py[Km_ind[i]] * part_py[Km_ind[i]] + part_pz[Km_ind[i]] * part_pz[Km_ind[i]]);
            p4_Km[i].SetPxPyPzE(part_px[Km_ind[i]], part_py[Km_ind[i]], part_pz[Km_ind[i]], sqrt(p * p + m_Km * m_Km));
            Km_FTOF_sec[i] = part_FTOF_sector_layer2[Km_ind[i]];
            if (abs(part_status[Km_ind[i]]) >= 2000 && abs(part_status[Km_ind[i]]) < 4000)
                Km_detect[i] = 2;
            if (part_status[Km_ind[i]] >= 4000)
                Km_detect[i] = 3;
            Km_chi2pid[i] = part_chi2pid[Km_ind[i]];
            if (part_DC_sector[Km_ind[i]] == 1)
                Km_sector[i] = 1;
            if (part_DC_sector[Km_ind[i]] == 2)
                Km_sector[i] = 2;
            if (part_DC_sector[Km_ind[i]] == 3)
                Km_sector[i] = 3;
            if (part_DC_sector[Km_ind[i]] == 4)
                Km_sector[i] = 4;
            if (part_DC_sector[Km_ind[i]] == 5)
                Km_sector[i] = 5;
            if (part_DC_sector[Km_ind[i]] == 6)
                Km_sector[i] = 6;
            // Km_RICHbest[i] = part_RICH_best[Km_ind[i]];
        }
    }
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// photon selector
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void select_photon(int run)
{

    float mom = 0;
    int check = 0;

    for (Int_t i = 0; i < Npart; i++)
    {
        if (i < BUFFER)
        {

            /// ////////////////////////////////////////////////////////////
            /// Forward detector cuts:

            if (abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000)
            {

                FD_photid_default_PID_check[i] = phot_default_PID_cut(i);
                FD_photid_charge_check[i] = phot_charge_cut(i);
                FD_photid_beta_check[i] = phot_beta_cut(i, run);
                // FD_photid_EC_sampling_fraction_check[i] = phot_EC_sampling_fraction_cut(i);
                FD_photid_EC_outer_vs_EC_inner_check[i] = phot_EC_outer_vs_EC_inner_cut(i);
                FD_photid_EC_hit_position_fiducial_check[i] = phot_EC_hit_position_fiducial_cut(i);

                if (FD_photid_default_PID_check[i])
                    FD_photid_default_PID_pass += 1;
                if (FD_photid_charge_check[i])
                    FD_photid_charge_pass += 1;
                if (FD_photid_beta_check[i] && FD_photid_charge_check[i] && FD_photid_default_PID_check[i])
                    FD_photid_beta_pass += 1;
                if (FD_photid_EC_hit_position_fiducial_check[i] && FD_photid_charge_check[i] && FD_photid_default_PID_check[i])
                    FD_photid_EC_hit_position_fiducial_pass += 1;

                if (use_own_PID_photon == true && simulation == false)
                {
                    if (FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_beta_check[i] && FD_photid_EC_hit_position_fiducial_check[i])
                    {
                        FD_photid_all_pass += 1;
                        FD_photid_all_check[i] = true;
                    }
                }

                if (use_own_PID_photon == true && simulation == true)
                {
                    if (FD_photid_default_PID_check[i] && FD_photid_charge_check[i] && FD_photid_EC_hit_position_fiducial_check[i])
                    {
                        FD_photid_all_pass += 1;
                        FD_photid_all_check[i] = true;
                    }
                }

                if (use_own_PID_photon == false)
                {
                    if (FD_photid_default_PID_check[i])
                    {
                        FD_photid_all_pass += 1;
                        FD_photid_all_check[i] = true;
                    }
                }
            }

            /// /////////////////////////////////////////////////////////////
            /// Forward tagger cuts:

            if (part_status[i] >= 1000 && part_status[i] < 2000)
            {

                FT_photid_PID_check[i] = FT_photid_PID_cut(i);
                FT_photid_charge_check[i] = FT_photid_charge_cut(i);
                FT_photid_FTCAL_fiducial_check[i] = FT_photid_FTCAL_fiducial_cut(i);
                FT_photid_beta_check[i] = FT_photid_beta_cut(i, run);

                if (FT_photid_PID_check[i])
                    FT_photid_PID_pass += 1;
                if (FT_photid_charge_check[i])
                    FT_photid_charge_pass += 1;
                if (FT_photid_FTCAL_fiducial_check[i] && FT_photid_charge_check[i] && FT_photid_PID_check[i])
                    FT_photid_FTCAL_fiducial_pass += 1;
                if (FT_photid_beta_check[i] && FT_photid_charge_check[i] && FT_photid_PID_check[i])
                    FT_photid_beta_pass += 1;

                if (FT_photid_PID_check[i] && FT_photid_FTCAL_fiducial_check[i])
                {
                    FT_photid_all_pass += 1;
                    FT_photid_all_check[i] = true;
                }
            }

            /// ///////////////////////////////////////////////////////////////
            /// Create photon selector

            bool selector;

            /// ///////////////////////////////////////////////////////////////
            /// Pick photon particle index and sort by momentum

            if (i < BUFFER)
            {
                mom = 0;
                check = 0;
                for (int k = 0; k < Npart; k++)
                {
                    if (k < BUFFER)
                    {

                        selector = (FD_photid_all_check[k] && abs(part_status[i]) >= 2000 && abs(part_status[i]) < 4000) || (FT_photid_all_check[k] && abs(part_status[i]) >= 1000 && abs(part_status[i]) < 2000);

                        if (selector)
                        { // photons only properly detected in FT and FD
                            for (int j = 0; j < i; j++)
                            {
                                if (k == g_ind[j])
                                    check = -1;
                            }

                            if (sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]) > mom && check != -1)
                            {
                                mom = sqrt(part_px[k] * part_px[k] + part_py[k] * part_py[k] + part_pz[k] * part_pz[k]);
                                g_ind[i] = k;
                                g_count += 1;
                            }
                            check = 0;
                        }
                    }
                }
            }
        }
    }

    /// //////////////////////////////////////////////////////////////////
    /// Assign properties

    for (int i = 0; i < BUFFER; i++)
    {
        if (g_ind[i] != -1)
        {
            g_vx[i] = part_vx[g_ind[i]];
            g_vy[i] = part_vy[g_ind[i]];
            g_vz[i] = part_vz[g_ind[i]];
            double p = sqrt(part_px[g_ind[i]] * part_px[g_ind[i]] + part_py[g_ind[i]] * part_py[g_ind[i]] + part_pz[g_ind[i]] * part_pz[g_ind[i]]);
            p4_phot[i].SetPxPyPzE(part_px[g_ind[i]], part_py[g_ind[i]], part_pz[g_ind[i]], p);
            g_sec[i] = part_Cal_PCAL_sector[g_ind[i]];
            if (abs(part_status[g_ind[i]]) >= 1000 && abs(part_status[g_ind[i]]) < 2000)
                phot_detect[i] = 1;
            if (abs(part_status[g_ind[i]]) >= 2000 && abs(part_status[g_ind[i]]) < 4000)
                phot_detect[i] = 2;
            if (part_Cal_PCAL_sector[g_ind[i]] == 1)
                phot_sector[i] = 1;
            if (part_Cal_PCAL_sector[g_ind[i]] == 2)
                phot_sector[i] = 2;
            if (part_Cal_PCAL_sector[g_ind[i]] == 3)
                phot_sector[i] = 3;
            if (part_Cal_PCAL_sector[g_ind[i]] == 4)
                phot_sector[i] = 4;
            if (part_Cal_PCAL_sector[g_ind[i]] == 5)
                phot_sector[i] = 5;
            if (part_Cal_PCAL_sector[g_ind[i]] == 6)
                phot_sector[i] = 6;
        }
    }
}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// select generated particles:

void write_gen(int run)
{

    TLorentzVector p4_gen[BUFFER];
    int mother[BUFFER];
    double mom[BUFFER];

    for (Int_t i = 0; i < BUFFER; i++)
    {
        p4_gen[i].SetPxPyPzE(0, 0, 0, 0);
        mother[i] = 0;
        mom[i] = 0;
    }

    for (int i = 0; i < BUFFER; i++)
    {
      mom[i] = sqrt(partMC_px[i] * partMC_px[i] + partMC_py[i] * partMC_py[i] + partMC_pz[i] * partMC_pz[i]);
      mother[i] = partMC_pid[partMC_mother[i]-1];     
      
      if(partMC_pid[i] == 11 && mom[i] > 0 && mom[i] < Ebeam-0.1){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_e * m_e));
        p4_gen_ele_px.push_back(p4_gen[i].Px());
        p4_gen_ele_py.push_back(p4_gen[i].Py());
        p4_gen_ele_pz.push_back(p4_gen[i].Pz());
        p4_gen_ele_E.push_back(p4_gen[i].E());
        p4_gen_ele_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 2212 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_p * m_p));
        p4_gen_prot_px.push_back(p4_gen[i].Px());
        p4_gen_prot_py.push_back(p4_gen[i].Py());
        p4_gen_prot_pz.push_back(p4_gen[i].Pz());
        p4_gen_prot_E.push_back(p4_gen[i].E());
        p4_gen_prot_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 2112 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_n * m_n));
        p4_gen_neutr_px.push_back(p4_gen[i].Px());
        p4_gen_neutr_py.push_back(p4_gen[i].Py());
        p4_gen_neutr_pz.push_back(p4_gen[i].Pz());
        p4_gen_neutr_E.push_back(p4_gen[i].E());
        p4_gen_neutr_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 211 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_pip * m_pip));
        p4_gen_pip_px.push_back(p4_gen[i].Px());
        p4_gen_pip_py.push_back(p4_gen[i].Py());
        p4_gen_pip_pz.push_back(p4_gen[i].Pz());
        p4_gen_pip_E.push_back(p4_gen[i].E());
        p4_gen_pip_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == -211 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_pim * m_pim));
        p4_gen_pim_px.push_back(p4_gen[i].Px());
        p4_gen_pim_py.push_back(p4_gen[i].Py());
        p4_gen_pim_pz.push_back(p4_gen[i].Pz());
        p4_gen_pim_E.push_back(p4_gen[i].E());
        p4_gen_pim_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 321 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_Kp * m_Kp));
        p4_gen_Kp_px.push_back(p4_gen[i].Px());
        p4_gen_Kp_py.push_back(p4_gen[i].Py());
        p4_gen_Kp_pz.push_back(p4_gen[i].Pz());
        p4_gen_Kp_E.push_back(p4_gen[i].E());
        p4_gen_Kp_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == -321 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + m_Km * m_Km));
        p4_gen_Km_px.push_back(p4_gen[i].Px());
        p4_gen_Km_py.push_back(p4_gen[i].Py());
        p4_gen_Km_pz.push_back(p4_gen[i].Pz());
        p4_gen_Km_E.push_back(p4_gen[i].E());
        p4_gen_Km_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 22 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i]));
        p4_gen_phot_px.push_back(p4_gen[i].Px());
        p4_gen_phot_py.push_back(p4_gen[i].Py());
        p4_gen_phot_pz.push_back(p4_gen[i].Pz());
        p4_gen_phot_E.push_back(p4_gen[i].E());
        p4_gen_phot_motherpid.push_back(mother[i]);
      } 
      if(partMC_pid[i] == 111 && mom[i] > 0){
        p4_gen[i].SetPxPyPzE(partMC_px[i], partMC_py[i], partMC_pz[i], sqrt(mom[i] * mom[i] + 0.134977 * 0.134977));
        p4_gen_pi0_px.push_back(p4_gen[i].Px());
        p4_gen_pi0_py.push_back(p4_gen[i].Py());
        p4_gen_pi0_pz.push_back(p4_gen[i].Pz());
        p4_gen_pi0_E.push_back(p4_gen[i].E());
        p4_gen_pi0_motherpid.push_back(mother[i]);
      }      
    }
}




/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// particle ID cuts for the FD:
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// hit in the FTOF is required for electrons and charged hadrons (basic cut)  --> For now: Selection of good paddels included
///

bool basic_FTOF_cut(int j)
{

    bool paddel_select = false;
    // if(part_FTOF_component_layer2[j] > 18 && part_FTOF_component_layer2[j] < 55) paddel_select = true;
    paddel_select = true; // use all paddels

    if ((part_FTOF_sector_layer1[j] != 0 || part_FTOF_sector_layer2[j] != 0 || part_FTOF_sector_layer3[j] != 0) && paddel_select)
        return true;
    else
        return false;
}

/// ///////////////////////////////////////////////////////////////////////////
/// FD electrons:

// electron has to have a hit in the FTOF, since this determines the start time of the event.

bool ele_default_PID_cut(int j)
{
    if (part_pid[j] == 11)
        return true;
    else
        return false;
}

bool ele_charge_cut(int j)
{
    if (part_charge[j] == -1)
        return true;
    else
        return false;
}

bool CC_nphe_cut(int j)
{

    double nphe_min = 2;

    if (part_CC_HTCC_nphe[j] > nphe_min)
        return true;
    else
        return false;
}

// EC cuts

bool EC_outer_vs_EC_inner_cut(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = false;
    bool loose = true;
    //////////////////////////

    double edep_min;

    if (loose == true)
    {
        edep_min = 0.06;
    }
    if (medium == true)
    {
        edep_min = 0.07;
    }
    if (tight == true)
    {
        edep_min = 0.09;
    }

    if (part_Cal_PCAL_energy[j] > edep_min)
        return true;
    else
        return false;
}



bool EC_sampling_fraction_cut(int j)
{

  // pass 2 total SF parameters (array entries represent sectors)
  
  //fall2018 inb:
  double p0mean_inb[] = {0.111767, 0.116619, 0.114606, 0.116586, 0.118251, 0.117391};
  double p1mean_inb[] = {-0.0281943, 0.0662751, -0.0896597, 0.181465, 0.085993, 0.0186504};
  double p2mean_inb[] = {0.00711137, 0.00633334, 0.00912098, 0.00652068, 0.00416682, 0.00622289};
  double p3mean_inb[] = {-0.000878776, -0.000780257, -0.00108891, -0.000645957, -0.000485189, -0.000829729};
  double p0sigma_inb[] = {-0.00497609, 0.0259435, 0.0296159, 0.0161445, 0.0239166, 0.0244309};
  double p1sigma_inb[] = {0.0275006, -0.000805156, -0.00449379, 0.0099462, 0.00192551, 0.00258059};
  double p2sigma_inb[] = {0.00253641, -0.00386759, -0.00469883, -0.00182968, -0.00355973, -0.00398967};
  double p3sigma_inb[] = {-0.000173549, 0.00030325, 0.000380195, 0.00012328, 0.000302528, 0.000340911};

  //fall2018 outb:
  double p0mean_outb[] = {0.111919, 0.11244, 0.11457, 0.124517, 0.109132, 0.115026};
  double p1mean_outb[] = {-0.00764253, 0.156704, 0.246338, 0.880436, -0.181137, 0.335205};
  double p2mean_outb[] = {0.00937217, 0.00924749, 0.00931085, 0.000420224, 0.00935877, 0.00756135};
  double p3mean_outb[] = {-0.000948603, -0.00095255, -0.000993028, -0.000123388, -0.000825096, -0.000822597};
  double p0sigma_outb[] = {-0.000828514, 0.019356, 0.023144, -0.000468566, 0.00500942, -0.00167471};
  double p1sigma_outb[] = {0.023998, 0.00249064, -0.00118722, 0.0223872, 0.0167932, 0.0243664};
  double p2sigma_outb[] = {0.000966279, -0.0022336, -0.00284786, 0.00143044, 0.000355477, 0.00124825};
  double p3sigma_outb[] = {-6.99914e-05, 0.000152666, 0.000164229, -0.000133009, -3.68797e-05, -0.000107538};

  //spring2019 inb:
  double p0mean_spr19[] = {0.11253, 0.113735, 0.112401, 0.115128, 0.113048, 0.1147};
  double p1mean_spr19[] = {-0.0689836, -0.044216, -0.160555, 0.108512, -0.153003, -0.0997027};
  double p2mean_spr19[] = {0.00793526, 0.00835112, 0.010781, 0.00695328, 0.00814596, 0.00766677};
  double p3mean_spr19[] = {-0.000920149, -0.000977549, -0.00128293, -0.000728472, -0.000957354, -0.00106035};
  double p0sigma_spr19[] = {0.0193473, 0.0351352, 0.0234448, 0.0238342, 0.0382829, 0.0125166};
  double p1sigma_spr19[] = {0.00436399, -0.0100767, 0.00133137, 0.00193097, -0.0127937, 0.0139726};
  double p2sigma_spr19[] = {-0.00222507, -0.00549204, -0.00353928, -0.00347335, -0.00623149, -0.00186593};
  double p3sigma_spr19[] = {0.000165603, 0.000440174, 0.000274543, 0.000260158, 0.000508396, 0.000206116};

  //MC inb:
  double p0mean_mcin[] = {0.118444, 0.118383, 0.118318, 0.118531, 0.117475, 0.119179};
  double p1mean_mcin[] = {-0.0445042, -0.0326496, 0.00402908, 0.0384926, -0.0768068, 0.0045002};
  double p2mean_mcin[] = {0.00522281, 0.00528503, 0.00542581, 0.00480286, 0.00601034, 0.00482612};
  double p3mean_mcin[] = {-0.000511299, -0.000504437, -0.000529496, -0.000449033, -0.000578208, -0.00046328};
  double p0sigma_mcin[] = {0.0204537, 0.0242836, 0.0320663, 0.0171258, 0.0236728, 0.0157762};
  double p1sigma_mcin[] = {0.00334198, -0.000192939, -0.00718171, 0.00698743, 0.000511942, 0.00801015};
  double p2sigma_mcin[] = {-0.00299872, -0.00369136, -0.00537603, -0.00246959, -0.00371171, -0.00215882};
  double p3sigma_mcin[] = {0.000232288, 0.000277151, 0.000425544, 0.00018372, 0.000290388, 0.000160679};

  //MC outb:
  double p0mean_mcout[] = {0.124075, 0.124086, 0.124071, 0.125947, 0.120091, 0.124457};
  double p1mean_mcout[] = {0.259169, 0.257774, 0.286948, 0.449492, -0.0180208, 0.274625};
  double p2mean_mcout[] = {0.00203391, 0.00204442, 0.00199255, 0.00074766, 0.00414372, 0.00178793};
  double p3mean_mcout[] = {-0.000135587, -0.000134399, -0.000132802, -3.89331e-05, -0.000283801, -0.000114508};
  double p0sigma_mcout[] = {0.00231891, 0.000686535, 0.000327404, -0.000165373, 0.00376051, 0.000856615};
  double p1sigma_mcout[] = {0.0205562, 0.0222621, 0.0226224, 0.0233338, 0.0190044, 0.0220325};
  double p2sigma_mcout[] = {0.000189274, 0.00044114, 0.000487635, 0.000531129, -9.69477e-06, 0.000416313};
  double p3sigma_mcout[] = {-3.11787e-05, -4.52271e-05, -4.89526e-05, -5.22239e-05, -2.02358e-05, -4.52287e-05}; 
    

  double sigma_range = 3.5;
  double mean = 0;
  double sigma = 0;
  double upper_lim_total = 0;
  double lower_lim_total = 0;
    
    for (Int_t k = 0; k < 6; k++)
    {
        if (part_Cal_PCAL_sector[j] - 1 == k)
        {
            if(inbending == true  && spring2019 == false && simulation == false){
                mean = p0mean_inb[k]*(1 + part_p[j]/sqrt(part_p[j]*part_p[j]+p1mean_inb[k])) + p2mean_inb[k]*part_p[j] + p3mean_inb[k]*part_p[j]*part_p[j];
                sigma = p0sigma_inb[k] + p1sigma_inb[k]/sqrt(part_p[j]) + p2sigma_inb[k]*part_p[j] + p3sigma_inb[k]*part_p[j]*part_p[j];
            }
            if(outbending == true  && simulation == false){
                mean = p0mean_outb[k]*(1 + part_p[j]/sqrt(part_p[j]*part_p[j]+p1mean_outb[k])) + p2mean_outb[k]*part_p[j] + p3mean_outb[k]*part_p[j]*part_p[j];
                sigma = p0sigma_outb[k] + p1sigma_outb[k]/sqrt(part_p[j]) + p2sigma_outb[k]*part_p[j] + p3sigma_outb[k]*part_p[j]*part_p[j];
            }
            if(spring2019 == true  && inbending == true && simulation == false){
                mean = p0mean_spr19[k]*(1 + part_p[j]/sqrt(part_p[j]*part_p[j]+p1mean_spr19[k])) + p2mean_spr19[k]*part_p[j] + p3mean_spr19[k]*part_p[j]*part_p[j];
                sigma = p0sigma_spr19[k] + p1sigma_spr19[k]/sqrt(part_p[j]) + p2sigma_spr19[k]*part_p[j] + p3sigma_spr19[k]*part_p[j]*part_p[j];
            } 
            if(inbending == true && simulation == true){
                mean = p0mean_mcin[k]*(1 + part_p[j]/sqrt(part_p[j]*part_p[j]+p1mean_mcin[k])) + p2mean_mcin[k]*part_p[j] + p3mean_mcin[k]*part_p[j]*part_p[j];
                sigma = p0sigma_mcin[k] + p1sigma_mcin[k]/sqrt(part_p[j]) + p2sigma_mcin[k]*part_p[j] + p3sigma_mcin[k]*part_p[j]*part_p[j];
            } 
            if(outbending == true  && simulation == true){
                mean = p0mean_mcout[k]*(1 + part_p[j]/sqrt(part_p[j]*part_p[j]+p1mean_mcout[k])) + p2mean_mcout[k]*part_p[j] + p3mean_mcout[k]*part_p[j]*part_p[j];
                sigma = p0sigma_mcout[k] + p1sigma_mcout[k]/sqrt(part_p[j]) + p2sigma_mcout[k]*part_p[j] + p3sigma_mcout[k]*part_p[j]*part_p[j];
            }     
            upper_lim_total = mean + sigma_range * sigma;
            lower_lim_total = mean - sigma_range * sigma;
        }
    }
    
    bool pass_band = part_Cal_energy_total[j] / part_p[j] <= upper_lim_total && part_Cal_energy_total[j] / part_p[j] >= lower_lim_total;
    
    /// ///
    /// triangle cut on SF PCAL vs SF ECin (array entries are momentum bins):

    // fall2018 inb:
    double p0_sec1_inb[] = {1.41582, 1.39934, 1.41204, 1.46385, 1.55892, 1.55892, 1.55892, 1.55892};
    double p1_sec1_inb[] = {0.212225, 0.215542, 0.217, 0.218279, 0.219881, 0.219881, 0.219881, 0.219881};
    double p0_sec2_inb[] = {1.44726, 1.44245, 1.47269, 1.53225, 1.61465, 1.61465, 1.61465, 1.61465};
    double p1_sec2_inb[] = {0.221991, 0.225772, 0.227888, 0.229099, 0.228898, 0.228898, 0.228898, 0.228898};
    double p0_sec3_inb[] = {1.38589, 1.3908, 1.42501, 1.48177, 1.57636, 1.57636, 1.57636, 1.57636};
    double p1_sec3_inb[] = {0.221492, 0.225738, 0.227955, 0.228604, 0.22836, 0.22836, 0.22836, 0.22836};
    double p0_sec4_inb[] = {1.38631, 1.38107, 1.39757, 1.44579, 1.54154, 1.54154, 1.54154, 1.54154};
    double p1_sec4_inb[] = {0.215784, 0.221511, 0.224982, 0.227812, 0.231076, 0.231076, 0.231076, 0.231076};
    double p0_sec5_inb[] = {1.50251, 1.52408, 1.52996, 1.49583, 1.39339, 1.39339, 1.39339, 1.39339};
    double p1_sec5_inb[] = {0.22202, 0.227163, 0.228794, 0.226487, 0.218168, 0.218168, 0.218168, 0.218168};
    double p0_sec6_inb[] = {1.51312, 1.52784, 1.57519, 1.67332, 1.85128, 1.85128, 1.85128, 1.85128};
    double p1_sec6_inb[] = {0.223651, 0.228082, 0.2305, 0.23241, 0.234238, 0.234238, 0.234238, 0.234238};

    // fall2018 outb:
    double p0_sec1_outb[] = {1.35967, 1.33697, 1.34111, 1.39563, 1.49066, 1.49066, 1.49066, 1.49066};
    double p1_sec1_outb[] = {0.21934, 0.222755, 0.224377, 0.227803, 0.23137, 0.23137, 0.23137, 0.23137};
    double p0_sec2_outb[] = {1.36974, 1.36895, 1.39344, 1.46945, 1.61251, 1.61251, 1.61251, 1.61251};
    double p1_sec2_outb[] = {0.218443, 0.223847, 0.227189, 0.231625, 0.238017, 0.238017, 0.238017, 0.238017};
    double p0_sec3_outb[] = {1.31891, 1.30602, 1.33372, 1.41351, 1.54453, 1.54453, 1.54453, 1.54453};
    double p1_sec3_outb[] = {0.219256, 0.223559, 0.22568, 0.230189, 0.235643, 0.235643, 0.235643, 0.235643};
    double p0_sec4_outb[] = {1.34425, 1.31016, 1.28753, 1.29671, 1.33831, 1.33831, 1.33831, 1.33831};
    double p1_sec4_outb[] = {0.217914, 0.221843, 0.221103, 0.221309, 0.222932, 0.222932, 0.222932, 0.222932};
    double p0_sec5_outb[] = {1.42433, 1.42395, 1.40932, 1.40816, 1.39868, 1.39868, 1.39868, 1.39868};
    double p1_sec5_outb[] = {0.218131, 0.222749, 0.225099, 0.225572, 0.224091, 0.224091, 0.224091, 0.224091};
    double p0_sec6_outb[] = {1.43741, 1.41924, 1.43218, 1.51807, 1.64554, 1.64554, 1.64554, 1.64554};
    double p1_sec6_outb[] = {0.220976, 0.225786, 0.228382, 0.232594, 0.238174, 0.238174, 0.238174, 0.238174};

    // spring 2019:
    double p0_sec1_spr19[] = {1.39979, 1.38131, 1.40144, 1.45945, 1.57144, 1.57144, 1.57144, 1.57144};
    double p1_sec1_spr19[] = {0.21633, 0.219878, 0.222027, 0.223849, 0.22592, 0.22592, 0.22592, 0.22592};
    double p0_sec2_spr19[] = {1.55243, 1.55871, 1.61819, 1.72756, 1.89424, 1.89424, 1.89424, 1.89424};
    double p1_sec2_spr19[] = {0.223711, 0.228475, 0.231997, 0.235231, 0.236813, 0.236813, 0.236813, 0.236813};
    double p0_sec3_spr19[] = {1.36528, 1.37519, 1.42453, 1.50339, 1.60863, 1.60863, 1.60863, 1.60863};
    double p1_sec3_spr19[] = {0.22, 0.224385, 0.227207, 0.228813, 0.228142, 0.228142, 0.228142, 0.228142};
    double p0_sec4_spr19[] = {1.38535, 1.3697, 1.39661, 1.4662, 1.63342, 1.63342, 1.63342, 1.63342};
    double p1_sec4_spr19[] = {0.214364, 0.219152, 0.222319, 0.225052, 0.229045, 0.229045, 0.229045, 0.229045};
    double p0_sec5_spr19[] = {1.47796, 1.47884, 1.48836, 1.47034, 1.38339, 1.38339, 1.38339, 1.38339};
    double p1_sec5_spr19[] = {0.219635, 0.223343, 0.224231, 0.221401, 0.210005, 0.210005, 0.210005, 0.210005};
    double p0_sec6_spr19[] = {1.51755, 1.53504, 1.59927, 1.73397, 2.00518, 2.00518, 2.00518, 2.00518};
    double p1_sec6_spr19[] = {0.222143, 0.226195, 0.228638, 0.23086, 0.233709, 0.233709, 0.233709, 0.233709};

    // mc inb:
    double p0_sec1_mcin[] = {1.30263, 1.30977, 1.31412, 1.31338, 1.31648, 1.31648, 1.31648, 1.31648};
    double p1_sec1_mcin[] = {0.224593, 0.227935, 0.229518, 0.22978, 0.22924, 0.22924, 0.22924, 0.22924};
    double p0_sec2_mcin[] = {1.30425, 1.30269, 1.30435, 1.30472, 1.29839, 1.29839, 1.29839, 1.29839};
    double p1_sec2_mcin[] = {0.224099, 0.226717, 0.22775, 0.227332, 0.224101, 0.224101, 0.224101, 0.224101};
    double p0_sec3_mcin[] = {1.28919, 1.28002, 1.24232, 1.18734, 1.12886, 1.12886, 1.12886, 1.12886};
    double p1_sec3_mcin[] = {0.223284, 0.225601, 0.224208, 0.220147, 0.214082, 0.214082, 0.214082, 0.214082};
    double p0_sec4_mcin[] = {1.29025, 1.29296, 1.29721, 1.30699, 1.32036, 1.32036, 1.32036, 1.32036};
    double p1_sec4_mcin[] = {0.221678, 0.224968, 0.226676, 0.227668, 0.227833, 0.227833, 0.227833, 0.227833};
    double p0_sec5_mcin[] = {1.2883, 1.27655, 1.24203, 1.19403, 1.15325, 1.15325, 1.15325, 1.15325};
    double p1_sec5_mcin[] = {0.224075, 0.226291, 0.225091, 0.22173, 0.217695, 0.217695, 0.217695, 0.217695};
    double p0_sec6_mcin[] = {1.30049, 1.30766, 1.31197, 1.31372, 1.31596, 1.31596, 1.31596, 1.31596};
    double p1_sec6_mcin[] = {0.224591, 0.227901, 0.229503, 0.23, 0.229422, 0.229422, 0.229422, 0.229422};

    // mc outb:

    double p0_sec1_mcout[] = {1.34088, 1.34685, 1.34914, 1.3544, 1.35628, 1.35628, 1.35628, 1.35628};
    double p1_sec1_mcout[] = {0.229433, 0.232865, 0.234838, 0.236458, 0.237532, 0.237532, 0.237532, 0.237532};
    double p0_sec2_mcout[] = {1.33865, 1.34439, 1.34938, 1.35247, 1.35423, 1.35423, 1.35423, 1.35423};
    double p1_sec2_mcout[] = {0.229232, 0.232698, 0.234938, 0.236421, 0.237479, 0.237479, 0.237479, 0.237479};
    double p0_sec3_mcout[] = {1.33981, 1.3419, 1.3459, 1.3506, 1.35108, 1.35108, 1.35108, 1.35108};
    double p1_sec3_mcout[] = {0.228797, 0.232055, 0.234322, 0.235882, 0.236784, 0.236784, 0.236784, 0.236784};
    double p0_sec4_mcout[] = {1.34856, 1.34544, 1.34355, 1.34317, 1.34383, 1.34383, 1.34383, 1.34383};
    double p1_sec4_mcout[] = {0.22849, 0.231908, 0.233365, 0.234226, 0.235005, 0.235005, 0.235005, 0.235005};
    double p0_sec5_mcout[] = {1.33882, 1.34416, 1.34414, 1.34681, 1.34737, 1.34737, 1.34737, 1.34737};
    double p1_sec5_mcout[] = {0.228946, 0.232535, 0.234382, 0.235812, 0.236718, 0.236718, 0.236718, 0.236718};
    double p0_sec6_mcout[] = {1.34048, 1.3457, 1.34985, 1.35331, 1.35303, 1.35303, 1.35303, 1.35303};
    double p1_sec6_mcout[] = {0.229432, 0.232767, 0.234884, 0.236384, 0.237264, 0.237264, 0.237264, 0.237264};

    int p_par = 0;
    if(part_p[j] <= 3)                  p_par = 0;  
    if(part_p[j] > 3 && part_p[j] <= 4) p_par = 1;
    if(part_p[j] > 4 && part_p[j] <= 5) p_par = 2; 
    if(part_p[j] > 5 && part_p[j] <= 6) p_par = 3; 
    if(part_p[j] > 6 && part_p[j] <= 7) p_par = 4; 
    if(part_p[j] > 7 && part_p[j] <= 8) p_par = 5; 
    if(part_p[j] > 8 && part_p[j] <= 9) p_par = 6; 
    if(part_p[j] > 9 )                  p_par = 7;  
  
    bool pass_triangle = false;
  
    if(inbending == true  && spring2019 == false && simulation == false){
      if(part_Cal_PCAL_sector[j] == 1) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec1_inb[p_par] - p0_sec1_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 2) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec2_inb[p_par] - p0_sec2_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 3) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec3_inb[p_par] - p0_sec3_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 4) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec4_inb[p_par] - p0_sec4_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 5) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec5_inb[p_par] - p0_sec5_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 6) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec6_inb[p_par] - p0_sec6_inb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
    }
    if(outbending == true  && simulation == false){
      if(part_Cal_PCAL_sector[j] == 1) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec1_outb[p_par] - p0_sec1_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 2) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec2_outb[p_par] - p0_sec2_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 3) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec3_outb[p_par] - p0_sec3_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 4) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec4_outb[p_par] - p0_sec4_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 5) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec5_outb[p_par] - p0_sec5_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 6) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec6_outb[p_par] - p0_sec6_outb[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
    }
    if(spring2019 == true  && inbending == true && simulation == false){
      if(part_Cal_PCAL_sector[j] == 1) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec1_spr19[p_par] - p0_sec1_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 2) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec2_spr19[p_par] - p0_sec2_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 3) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec3_spr19[p_par] - p0_sec3_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 4) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec4_spr19[p_par] - p0_sec4_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 5) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec5_spr19[p_par] - p0_sec5_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 6) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec6_spr19[p_par] - p0_sec6_spr19[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
    } 
    if(inbending == true && simulation == true){
      if(part_Cal_PCAL_sector[j] == 1) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec1_mcin[p_par] - p0_sec1_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 2) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec2_mcin[p_par] - p0_sec2_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 3) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec3_mcin[p_par] - p0_sec3_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 4) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec4_mcin[p_par] - p0_sec4_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 5) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec5_mcin[p_par] - p0_sec5_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 6) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec6_mcin[p_par] - p0_sec6_mcin[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
    } 
    if(outbending == true  && simulation == true){
      if(part_Cal_PCAL_sector[j] == 1) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec1_mcout[p_par] - p0_sec1_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 2) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec2_mcout[p_par] - p0_sec2_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 3) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec3_mcout[p_par] - p0_sec3_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 4) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec4_mcout[p_par] - p0_sec4_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 5) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec5_mcout[p_par] - p0_sec5_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
      if(part_Cal_PCAL_sector[j] == 6) pass_triangle = part_Cal_PCAL_energy[j] / part_p[j] > (p1_sec6_mcout[p_par] - p0_sec6_mcout[p_par]*part_Cal_ECin_energy[j] / part_p[j]);
    }     

    bool pass_threshold = (part_Cal_PCAL_energy[j] / part_p[j]) > 0.05;

    // old simple pass1 cut:
    //if(part_p[j] < 4.5){pass_triangle = true;}
    //else{ pass_triangle = part_Cal_ECin_energy[j] / part_p[j] > (0.2 - part_Cal_PCAL_energy[j] / part_p[j]);}
    
    ///
    /// ///

    if(pass_band && pass_triangle && pass_threshold) return true;
    else return false;
}





/////////////////////////////////////////////////////////////
/// EC hit position homogenous cut

bool EC_hit_position_fiducial_cut_homogeneous(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = false;
    bool loose = true;
    //////////////////////////

    // Cut using the natural directions of the scintillator bars/ fibers:

    double v = part_Cal_PCAL_lv[j];
    double w = part_Cal_PCAL_lw[j];

    /// v + w is going from the side to the back end of the PCAL, u is going from side to side
    /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.
    /// a cut is only applied on v and w

    ///////////////////////////////////////////////////////////////////
    /// inbending:
    //
    double min_v_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_v_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_v_loose_inb[] = {9.0, 9.0, 9.0, 13.5, 9.0, 9.0};
    //
    double max_v_tight_inb[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_inb[] = {400, 400, 400, 400, 400, 400};
    double max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
    //
    double min_w_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_w_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_w_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_w_tight_inb[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_inb[] = {400, 400, 400, 400, 400, 400};
    double max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};

    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    //
    double min_v_tight_out[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_v_med_out[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_v_loose_out[] = {9.0, 9.0, 9.0, 13.5, 9.0, 9.0};
    //
    double max_v_tight_out[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_out[] = {400, 400, 400, 400, 400, 400};
    double max_v_loose_out[] = {400, 400, 400, 400, 400, 400};
    //
    double min_w_tight_out[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_w_med_out[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_w_loose_out[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_w_tight_out[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_out[] = {400, 400, 400, 400, 400, 400};
    double max_w_loose_out[] = {400, 400, 400, 400, 400, 400};

    //////////////////////////////////////////////////////////////

    double min_v = 0;
    double max_v = 0;
    double min_w = 0;
    double max_w = 0;

    for (Int_t k = 0; k < 6; k++)
    {
        if (part_Cal_PCAL_sector[j] - 1 == k && inbending == true)
        {
            if (tight == true)
            {
                min_v = min_v_tight_inb[k];
                max_v = max_v_tight_inb[k];
                min_w = min_w_tight_inb[k];
                max_w = max_w_tight_inb[k];
            }
            if (medium == true)
            {
                min_v = min_v_med_inb[k];
                max_v = max_v_med_inb[k];
                min_w = min_w_med_inb[k];
                max_w = max_w_med_inb[k];
            }
            if (loose == true)
            {
                min_v = min_v_loose_inb[k];
                max_v = max_v_loose_inb[k];
                min_w = min_w_loose_inb[k];
                max_w = max_w_loose_inb[k];
            }
        }
        if (part_Cal_PCAL_sector[j] - 1 == k && outbending == true)
        {
            if (tight == true)
            {
                min_v = min_v_tight_out[k];
                max_v = max_v_tight_out[k];
                min_w = min_w_tight_out[k];
                max_w = max_w_tight_out[k];
            }
            if (medium == true)
            {
                min_v = min_v_med_out[k];
                max_v = max_v_med_out[k];
                min_w = min_w_med_out[k];
                max_w = max_w_med_out[k];
            }
            if (loose == true)
            {
                min_v = min_v_loose_out[k];
                max_v = max_v_loose_out[k];
                min_w = min_w_loose_out[k];
                max_w = max_w_loose_out[k];
            }
        }
    }

    if (v > min_v && v < max_v && w > min_w && w < max_w)
        return true;
    else
        return false;
}



/// //////////////////////////////////////////////////////////////////
/// edge based RG-A DC fiducial cuts (pass2):


bool DC_fiducial_cut_edge(int j, int region)
{

  double DCedge_ele_inb[]  = {5.0, 5.0, 10.0};
  double DCedge_prot_inb[] = {2.5, 2.5, 9.0};
  double DCedge_pip_inb[]  = {2.5, 2.5, 9.0};
  double DCedge_pim_inb[]  = {3.5, 3.0, 7.0};
  double DCedge_Kp_inb[]   = {2.5, 2.0, 9.0};
  double DCedge_Km_inb[]   = {3.5, 2.5, 5.0};

  double DCedge_ele_outb[]  = {3.0, 3.0, 10.0};
  double DCedge_prot_outb[] = {3.5, 3.0, 7.0};
  double DCedge_pip_outb[]  = {3.5, 2.5, 6.5};
  double DCedge_pim_outb[]  = {2.5, 2.5, 10.0};
  double DCedge_Kp_outb[]   = {3.5, 2.5, 6.5};
  double DCedge_Km_outb[]   = {2.5, 2.5, 10.0};


  double edge_cut = 0;

  if(outbending == false && region == 1 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_inb[2];
  if(outbending == false && region == 1 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_inb[2];
  if(outbending == false && region == 1 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_inb[2];
  if(outbending == false && region == 1 &&  part_pid[j] == -211) edge_cut = DCedge_pim_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == -211) edge_cut = DCedge_pim_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == -211) edge_cut = DCedge_pim_inb[2];
  if(outbending == false && region == 1 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_inb[2];
  if(outbending == false && region == 1 &&  part_pid[j] == -321) edge_cut = DCedge_Km_inb[0];
  if(outbending == false && region == 2 &&  part_pid[j] == -321) edge_cut = DCedge_Km_inb[1];
  if(outbending == false && region == 3 &&  part_pid[j] == -321) edge_cut = DCedge_Km_inb[2];
  
  if(outbending == true && region == 1 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == 11)   edge_cut = DCedge_ele_outb[2];
  if(outbending == true && region == 1 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == 2212) edge_cut = DCedge_prot_outb[2];
  if(outbending == true && region == 1 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == 211)  edge_cut = DCedge_pip_outb[2];
  if(outbending == true && region == 1 &&  part_pid[j] == -211) edge_cut = DCedge_pim_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == -211) edge_cut = DCedge_pim_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == -211) edge_cut = DCedge_pim_outb[2];
  if(outbending == true && region == 1 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == 321)  edge_cut = DCedge_Kp_outb[2];
  if(outbending == true && region == 1 &&  part_pid[j] == -321) edge_cut = DCedge_Km_outb[0];
  if(outbending == true && region == 2 &&  part_pid[j] == -321) edge_cut = DCedge_Km_outb[1];
  if(outbending == true && region == 3 &&  part_pid[j] == -321) edge_cut = DCedge_Km_outb[2];
  
  double edge_val;

  switch(region)
  {
    case 1:
      edge_val = part_DC_edge1[j];
      break;
    case 2:
      edge_val = part_DC_edge2[j];
      break;
    case 3:
      edge_val = part_DC_edge3[j];
      break;
    default:
      edge_val = 0;
      break;
  }
    
  if(edge_val > edge_cut) return true;
  else return false;
}


/////////////////////////////////////////////////////////////////////

bool DC_z_vertex_cut(int j)
{

    // pass 2 (adjusted to cross sections)

    double vz_min_sect_inb[] = {-8, -8, -8, -8, -8, -8};
    double vz_max_sect_inb[] = {2, 2, 2, 2, 2, 2};

    double vz_min_sect_outb[] = {-11, -11, -11, -11, -11, -11};
    double vz_max_sect_outb[] = {1, 1, 1, 1, 1, 1};

    double vz_min_sect[6];
    double vz_max_sect[6];

    for (Int_t i = 0; i < 6; i++)
    {
        if (inbending == true)
        {
            vz_min_sect[i] = vz_min_sect_inb[i];
            vz_max_sect[i] = vz_max_sect_inb[i];
        }
        if (outbending == true)
        {
            vz_min_sect[i] = vz_min_sect_outb[i];
            vz_max_sect[i] = vz_max_sect_outb[i];
        }
    }

    double vz_min = 0;
    double vz_max = 0;

    for (Int_t k = 0; k < 6; k++)
    {
        if (part_Cal_PCAL_sector[j] - 1 == k)
        {
            vz_min = vz_min_sect[k];
            vz_max = vz_max_sect[k];
        }
    }

    if (part_vz[j] > vz_min && part_vz[j] < vz_max)
        return true;
    else
        return false;
}



/// ///////////////////////////////////////////////////////////////////
/// FD hadrons:

// a) default:

bool prot_default_PID_cut(int j)
{
    if (part_pid[j] == 2212)
        return true;
    else
        return false;
}
bool neutr_default_PID_cut(int j)
{
    if (part_pid[j] == 2112)
        return true;
    else
        return false;
}
bool pip_default_PID_cut(int j)
{
    if (part_pid[j] == 211)
        return true;
    else
        return false;
}
bool pim_default_PID_cut(int j)
{
    if (part_pid[j] == -211)
        return true;
    else
        return false;
}
bool Kp_default_PID_cut(int j)
{
    if (part_pid[j] == 321)
        return true;
    else
        return false;
}
bool Km_default_PID_cut(int j)
{
    if (part_pid[j] == -321)
        return true;
    else
        return false;
}

// b) charge cuts

bool prot_charge_cut(int j)
{
    if (part_charge[j] == +1)
        return true;
    else
        return false;
}
bool neutr_charge_cut(int j)
{
    if (part_charge[j] == 0)
        return true;
    else
        return false;
}
bool pip_charge_cut(int j)
{
    if (part_charge[j] == +1)
        return true;
    else
        return false;
}
bool pim_charge_cut(int j)
{
    if (part_charge[j] == -1)
        return true;
    else
        return false;
}
bool Kp_charge_cut(int j)
{
    if (part_charge[j] == +1)
        return true;
    else
        return false;
}
bool Km_charge_cut(int j)
{
    if (part_charge[j] == -1)
        return true;
    else
        return false;
}

// b.1) electron rejection cut for negative hadrons

bool pim_ele_reject_cut(int j)
{
    if (FD_eid_all_check[j] == false)
        return true;
    else
        return false;
}
bool Km_ele_reject_cut(int j)
{
    if (FD_eid_all_check[j] == false)
        return true;
    else
        return false;
}

bool pim_EC_outer_vs_EC_inner_cut(int j)
{
    double edep_max = 0.06;
    if (part_Cal_PCAL_energy[j] < edep_max)
        return true;
    else
        return false;
}

bool Km_EC_outer_vs_EC_inner_cut(int j)
{
    double edep_max = 0.06;
    if (part_Cal_PCAL_energy[j] < edep_max)
        return true;
    else
        return false;
}

// DC cuts

bool DC_hit_position_region1_fiducial_cut_hadrons_positive(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 1.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -1.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {15, 15, 15, 15, 15, 15};
    double height_outb[] = {20, 20, 20, 20, 20, 20};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 25 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 30 + add;
        }
    }

    double x1_rot = part_DC_c1y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c1x[j] * cos(sec * 60.0 * Pival / 180);
    double y1_rot = part_DC_c1y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c1x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y1_rot);
    double right = (height + slope * y1_rot);

    double radius2_DCr1_min = pow(radius, 2) - pow(y1_rot, 2); // cut out the inner circle
    double radius2_DCr1_max = pow(155, 2) - pow(y1_rot, 2);    // cut out the outer circle

    if (x1_rot > left && x1_rot > right && pow(x1_rot, 2) > radius2_DCr1_min && pow(x1_rot, 2) < radius2_DCr1_max)
        return true;
    else
        return false;
}

bool DC_hit_position_region2_fiducial_cut_hadrons_positive(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 2.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -2.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {25, 25, 25, 25, 25, 25};
    double height_outb[] = {31, 31, 31, 31, 31, 31};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 40 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 50 + add;
        }
    }

    double x2_rot = part_DC_c2y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c2x[j] * cos(sec * 60.0 * Pival / 180);
    double y2_rot = part_DC_c2y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c2x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y2_rot);
    double right = (height + slope * y2_rot);

    double radius2_DCr2_min = pow(radius, 2) - pow(y2_rot, 2); // cut out the inner circle
    double radius2_DCr2_max = pow(245, 2) - pow(y2_rot, 2);    // cut out the outer circle

    // return true;
    if (x2_rot > left && x2_rot > right && pow(x2_rot, 2) > radius2_DCr2_min && pow(x2_rot, 2) < radius2_DCr2_max)
        return true;
    else
        return false;
}

bool DC_hit_position_region3_fiducial_cut_hadrons_positive(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 3.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -3.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {38, 38, 38, 38, 38, 38};
    double height_outb[] = {43, 43, 43, 43, 43, 43};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 52 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 62 + add;
        }
    }

    double x3_rot = part_DC_c3y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c3x[j] * cos(sec * 60.0 * Pival / 180);
    double y3_rot = part_DC_c3y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c3x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y3_rot);
    double right = (height + slope * y3_rot);

    double radius2_DCr3_min = pow(radius, 2) - pow(y3_rot, 2); // cut out the inner circle
    double radius2_DCr3_max = pow(355, 2) - pow(y3_rot, 2);    // cut out the outer circle

    // return true;
    if (x3_rot > left && x3_rot > right && pow(x3_rot, 2) > radius2_DCr3_min && pow(x3_rot, 2) < radius2_DCr3_max)
        return true;
    else
        return false;
}

bool DC_hit_position_region1_fiducial_cut_hadrons_negative(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 1.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -1.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {15, 15, 15, 15, 15, 15};
    double height_outb[] = {15, 15, 15, 15, 15, 15};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 24 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 25 + add;
        }
    }

    double x1_rot = part_DC_c1y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c1x[j] * cos(sec * 60.0 * Pival / 180);
    double y1_rot = part_DC_c1y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c1x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y1_rot);
    double right = (height + slope * y1_rot);

    double radius2_DCr1_min = pow(radius, 2) - pow(y1_rot, 2); // cut out the inner circle
    double radius2_DCr1_max = pow(155, 2) - pow(y1_rot, 2);    // cut out the outer circle

    if (x1_rot > left && x1_rot > right && pow(x1_rot, 2) > radius2_DCr1_min && pow(x1_rot, 2) < radius2_DCr1_max)
        return true;
    else
        return false;
}

bool DC_hit_position_region2_fiducial_cut_hadrons_negative(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 2.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -2.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {26, 26, 26, 26, 26, 26};
    double height_outb[] = {29, 29, 29, 29, 29, 29};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 42 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 43 + add;
        }
    }

    double x2_rot = part_DC_c2y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c2x[j] * cos(sec * 60.0 * Pival / 180);
    double y2_rot = part_DC_c2y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c2x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y2_rot);
    double right = (height + slope * y2_rot);

    double radius2_DCr2_min = pow(radius, 2) - pow(y2_rot, 2); // cut out the inner circle
    double radius2_DCr2_max = pow(245, 2) - pow(y2_rot, 2);    // cut out the outer circle

    // return true;
    if (x2_rot > left && x2_rot > right && pow(x2_rot, 2) > radius2_DCr2_min && pow(x2_rot, 2) < radius2_DCr2_max)
        return true;
    else
        return false;
}

bool DC_hit_position_region3_fiducial_cut_hadrons_negative(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double add = 0; // value in cm added to the height and radius of the cut
    if (tight == true)
    {
        add = 3.0;
    }
    if (medium == true)
    {
        add = 0.0;
    }
    if (loose == true)
    {
        add = -3.0;
    }

    double angle = 60;
    double height = 0;
    double radius = 0;

    double height_inb[] = {32, 32, 32, 32, 32, 32};
    double height_outb[] = {53, 53, 53, 53, 53, 53};

    int sec = part_DC_sector[j] - 1;

    for (Int_t k = 0; k < 6; k++)
    {
        if (sec == k && inbending == true)
        {
            height = height_inb[k] + add;
            radius = 46 + add;
        }
        if (sec == k && outbending == true)
        {
            height = height_outb[k] + add;
            radius = 68 + add;
        }
    }

    double x3_rot = part_DC_c3y[j] * sin(sec * 60.0 * Pival / 180) + part_DC_c3x[j] * cos(sec * 60.0 * Pival / 180);
    double y3_rot = part_DC_c3y[j] * cos(sec * 60.0 * Pival / 180) - part_DC_c3x[j] * sin(sec * 60.0 * Pival / 180);

    double slope = 1 / tan(0.5 * angle * Pival / 180);
    double left = (height - slope * y3_rot);
    double right = (height + slope * y3_rot);

    double radius2_DCr3_min = pow(radius, 2) - pow(y3_rot, 2); // cut out the inner circle
    double radius2_DCr3_max = pow(355, 2) - pow(y3_rot, 2);    // cut out the outer circle

    // return true;
    if (x3_rot > left && x3_rot > right && pow(x3_rot, 2) > radius2_DCr3_min && pow(x3_rot, 2) < radius2_DCr3_max)
        return true;
    else
        return false;
}



// g) delta vz cuts

bool prot_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool neutr_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool pip_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool pim_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool Kp_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool Km_delta_vz_cut(int j)
{

    double dvz_min = -20;
    double dvz_max = 20;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

// /////////////////////////////////////////////////////////////////
// FD photons:

bool phot_default_PID_cut(int j)
{
    if (part_pid[j] == 22)
        return true;
    else
        return false;
}

bool phot_charge_cut(int j)
{
    if (part_charge[j] == 0)
        return true;
    else
        return false;
}

bool phot_beta_cut(int j, int run)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    double min, max;

    if (loose == true)
    { // no additional cut
        min = 0.9;
        max = 2.0;
    }

    if (medium == true)
    { 
        min = 0.9;
        max = 1.1;
    }

    if (tight == true)
    { 
        min = 0.95;
        max = 1.05;
    }

    if (Beta_neutral(j, run) > min && Beta_neutral(j, run) < max && part_p[j] > 0.10)
        return true;

    else
        return false;
}

bool phot_EC_sampling_fraction_cut(int j)
{

    double ECfrac_min = 0;
    double ECfrac_max = 1;

    if (part_Cal_energy_total[j] / part_p[j] > ECfrac_min && part_Cal_energy_total[j] / part_p[j] < ECfrac_max)
        return true;
    else
        return false;
}

bool phot_EC_outer_vs_EC_inner_cut(int j)
{

    double edep_min = 0.01;

    if (part_Cal_ECin_energy[j] + part_Cal_ECout_energy[j] > edep_min)
        return true;
    else
        return false;
}

bool phot_EC_hit_position_fiducial_cut(int j)
{

    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    //////////////////////////

    // Cut using the natural directions of the scintillator bars/ fibers:

    double v = part_Cal_PCAL_lv[j];
    double w = part_Cal_PCAL_lw[j];

    /// v + w is going from the side to the back end of the PCAL, u is going from side to side
    /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.
    /// a cut is only applied on v and w

    ///////////////////////////////////////////////////////////////////
    /// inbending:
    //
    double min_v_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_v_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_v_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_v_tight_inb[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_inb[] = {400, 400, 400, 400, 400, 400};
    double max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
    //
    double min_w_tight_inb[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_w_med_inb[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_w_loose_inb[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_w_tight_inb[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_inb[] = {400, 400, 400, 400, 400, 400};
    double max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};

    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    //
    double min_v_tight_out[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_v_med_out[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_v_loose_out[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_v_tight_out[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_out[] = {400, 400, 400, 400, 400, 400};
    double max_v_loose_out[] = {400, 400, 400, 400, 400, 400};
    //
    double min_w_tight_out[] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
    double min_w_med_out[] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
    double min_w_loose_out[] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    //
    double max_w_tight_out[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_out[] = {400, 400, 400, 400, 400, 400};
    double max_w_loose_out[] = {400, 400, 400, 400, 400, 400};

    //////////////////////////////////////////////////////////////

    double min_v = 0;
    double max_v = 0;
    double min_w = 0;
    double max_w = 0;

    for (Int_t k = 0; k < 6; k++)
    {
        if (part_Cal_PCAL_sector[j] - 1 == k && inbending == true)
        {
            if (tight == true)
            {
                min_v = min_v_tight_inb[k];
                max_v = max_v_tight_inb[k];
                min_w = min_w_tight_inb[k];
                max_w = max_w_tight_inb[k];
            }
            if (medium == true)
            {
                min_v = min_v_med_inb[k];
                max_v = max_v_med_inb[k];
                min_w = min_w_med_inb[k];
                max_w = max_w_med_inb[k];
            }
            if (loose == true)
            {
                min_v = min_v_loose_inb[k];
                max_v = max_v_loose_inb[k];
                min_w = min_w_loose_inb[k];
                max_w = max_w_loose_inb[k];
            }
        }
        if (part_Cal_PCAL_sector[j] - 1 == k && outbending == true)
        {
            if (tight == true)
            {
                min_v = min_v_tight_out[k];
                max_v = max_v_tight_out[k];
                min_w = min_w_tight_out[k];
                max_w = max_w_tight_out[k];
            }
            if (medium == true)
            {
                min_v = min_v_med_out[k];
                max_v = max_v_med_out[k];
                min_w = min_w_med_out[k];
                max_w = max_w_med_out[k];
            }
            if (loose == true)
            {
                min_v = min_v_loose_out[k];
                max_v = max_v_loose_out[k];
                min_w = min_w_loose_out[k];
                max_w = max_w_loose_out[k];
            }
        }
    }

    if (v > min_v && v < max_v && w > min_w && w < max_w)
        return true;
    else
        return false;
}



/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Particle ID cuts for the Forward Tagger:
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// ////////////////////////////////////////////////////////
/// FT electrons:

bool FT_eid_charge_cut(int j)
{
    if (part_charge[j] == -1)
        return true;
    else
        return false;
}

bool FT_eid_PID_cut(int j)
{
    if (part_pid[j] == 11)
        return true;
    else
        return false;
}

bool FT_eid_FTCAL_fiducial_cut(int j)
{

    double p = sqrt(part_px[j] * part_px[j] + part_py[j] * part_py[j] + part_pz[j] * part_pz[j]);
    double theta = acos(part_pz[j] / p);

    if (theta * 180 / Pival > 2.5 && theta * 180 / Pival < 4.5)
        return true;
    else
        return false;
}

bool FT_eid_FTTRK_fiducial_cut(int j)
{

    double p = sqrt(part_px[j] * part_px[j] + part_py[j] * part_py[j] + part_pz[j] * part_pz[j]);
    double theta = acos(part_pz[j] / p);

    if (theta * 180 / Pival > 2.5 && theta * 180 / Pival < 4.5)
        return true;
    else
        return false;
}

bool FT_eid_FTHODO_fiducial_cut(int j)
{

    double p = sqrt(part_px[j] * part_px[j] + part_py[j] * part_py[j] + part_pz[j] * part_pz[j]);
    double theta = acos(part_pz[j] / p);

    if (theta * 180 / Pival > 2.5 && theta * 180 / Pival < 4.5)
        return true;
    else
        return false;
}

bool FT_eid_energy_vs_radius_cut(int j)
{

    return true;
}

/// /////////////////////////////////////////////////
/// FT photons

bool FT_photid_charge_cut(int j)
{
    if (part_charge[j] == 0)
        return true;
    else
        return false;
}

bool FT_photid_PID_cut(int j)
{
    if (part_pid[j] == 22)
        return true;
    else
        return false;
}

bool FT_photid_FTCAL_fiducial_cut(int j)
{

    TVector3 V3ECalPos(part_FTCAL_x[j], part_FTCAL_y[j], 0);
    double clusX = part_FTCAL_x[j];
    double clusY = part_FTCAL_y[j];

    bool res = true && V3ECalPos.Mag() > 8 && V3ECalPos.Mag() < 15 && TMath::Power(clusX + 8.5, 2) + TMath::Power(clusY - 10, 2) > 1.5 * 1.5 && TMath::Power(clusX + 10, 2) + TMath::Power(clusY + 5, 2) > 1.5 * 1.5 && TMath::Power(clusX + 6, 2) + TMath::Power(clusY + 13.5, 2) > 2 * 2 && TMath::Power(clusX - 4, 2) + TMath::Power(clusY + 6.7, 2) > 1.5 * 1.5 && TMath::Power(clusX - 6, 2) + TMath::Power(clusY + 6, 2) > 1;

    return res;
}

bool FT_photid_beta_cut(int j, int run)
{

    return true;
}

/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Particle ID cuts for the Central Detector:
/// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// beta cuts


bool neutr_beta_cut(int j, int run){
  if(Beta_neutral(j, run) < 0.95 && Beta_neutral(j, run) > 0) return true;
  else return false;
}


bool CD_neutr_beta_cut(int j, int run)
{
    if (Beta_neutral(j, run) < 0.95 && Beta_neutral(j, run) > 0)
        return true;
    else
        return false;
}

// dvz cuts

bool CD_prot_delta_vz_cut(int j)
{

    double mean = 0.7486;
    double sigma = 3.237;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool CD_neutr_delta_vz_cut(int j)
{

    double mean = 2.254;
    double sigma = 2.693;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool CD_pip_delta_vz_cut(int j)
{

    double mean = -0.6183;
    double sigma = 3.684;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool CD_pim_delta_vz_cut(int j)
{

    double mean = -0.5485;
    double sigma = 3.677;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool CD_Kp_delta_vz_cut(int j)
{

    double mean = 1.658;
    double sigma = 2.52;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}

bool CD_Km_delta_vz_cut(int j)
{

    double mean = -1.161;
    double sigma = 2.691;
    double dvz_min = mean - 3 * sigma;
    double dvz_max = mean + 3 * sigma;

    dvz_min = -20.0;
    dvz_max = 20.0;

    if (Getdvz(j) > dvz_min && Getdvz(j) < dvz_max)
        return true;
    else
        return false;
}




/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// fill output tree variables (particles from all detrectors are combined):
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_output_vector_electron(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_ele[i].E() > 0 && ele_detect[i] == 2)
        {
            p4_ele_px.push_back(p4_ele[i].Px());
            p4_ele_py.push_back(p4_ele[i].Py());
            p4_ele_pz.push_back(p4_ele[i].Pz());
            p4_ele_E.push_back(p4_ele[i].E());
            ele_det.push_back(ele_detect[i]);
            ele_sec.push_back(ele_sector[i]);
            p4_ele_chi2pid.push_back(ele_chi2pid[i]);
            if(fill_coord_ele == true){
               p4_ele_dcx1.push_back(ele_dcx1[i]); 
               p4_ele_dcy1.push_back(ele_dcy1[i]); 
               p4_ele_dcx2.push_back(ele_dcx2[i]); 
               p4_ele_dcy2.push_back(ele_dcy2[i]); 
               p4_ele_dcx3.push_back(ele_dcx3[i]); 
               p4_ele_dcy3.push_back(ele_dcy3[i]); 
               p4_ele_pcalv.push_back(ele_pcalv[i]); 
               p4_ele_pcalw.push_back(ele_pcalw[i]);
            } 
            if(fill_edge_ele == true){
               p4_ele_dcedge1.push_back(ele_dcedge1[i]); 
               p4_ele_dcedge2.push_back(ele_dcedge2[i]); 
               p4_ele_dcedge3.push_back(ele_dcedge3[i]);  
            }      
        }
    }
}

void fill_output_vector_proton(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_prot[i].E() > 0)
        {
            p4_prot_px.push_back(p4_prot[i].Px());
            p4_prot_py.push_back(p4_prot[i].Py());
            p4_prot_pz.push_back(p4_prot[i].Pz());
            p4_prot_E.push_back(p4_prot[i].E());
            prot_det.push_back(prot_detect[i]);
            prot_sec.push_back(prot_sector[i]);
            p4_prot_chi2pid.push_back(prot_chi2pid[i]);
            if(fill_coord_prot == true){
              p4_prot_dcx1.push_back(prot_dcx1[i]); 
              p4_prot_dcy1.push_back(prot_dcy1[i]); 
              p4_prot_dcz1.push_back(prot_dcz1[i]); 
            }
            if(fill_edge_prot == true){
               p4_prot_dcedge1.push_back(prot_dcedge1[i]); 
               p4_prot_dcedge2.push_back(prot_dcedge2[i]); 
               p4_prot_dcedge3.push_back(prot_dcedge3[i]); 
            }  
        }
    }
}

void fill_output_vector_neutron(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_neutr[i].E() > 0)
        {
            p4_neutr_px.push_back(p4_neutr[i].Px());
            p4_neutr_py.push_back(p4_neutr[i].Py());
            p4_neutr_pz.push_back(p4_neutr[i].Pz());
            p4_neutr_E.push_back(p4_neutr[i].E());
            neutr_det.push_back(neutr_detect[i]);
            neutr_sec.push_back(neutr_sector[i]);
            p4_neutr_chi2pid.push_back(neutr_chi2pid[i]);
        }
    }
}

void fill_output_vector_pip(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_pip[i].E() > 0)
        {
            p4_pip_px.push_back(p4_pip[i].Px());
            p4_pip_py.push_back(p4_pip[i].Py());
            p4_pip_pz.push_back(p4_pip[i].Pz());
            p4_pip_E.push_back(p4_pip[i].E());
            pip_det.push_back(pip_detect[i]);
            pip_sec.push_back(pip_sector[i]);
            p4_pip_chi2pid.push_back(pip_chi2pid[i]);   
            if(fill_coord_pip == true){
              p4_pip_dcx1.push_back(pip_dcx1[i]); 
              p4_pip_dcy1.push_back(pip_dcy1[i]); 
              p4_pip_dcx2.push_back(pip_dcx2[i]); 
              p4_pip_dcy2.push_back(pip_dcy2[i]); 
              p4_pip_dcx3.push_back(pip_dcx3[i]); 
              p4_pip_dcy3.push_back(pip_dcy3[i]);  
            } 
            if(fill_edge_pip == true){
               p4_pip_dcedge1.push_back(pip_dcedge1[i]); 
               p4_pip_dcedge2.push_back(pip_dcedge2[i]); 
               p4_pip_dcedge3.push_back(pip_dcedge3[i]); 
            }  
            p4_pip_ml.push_back(ml_value(pip_ind[i]));
            p4_pip_beta.push_back(pip_beta[i]);
            p4_pip_LTCC_nphe.push_back(part_CC_LTCC_nphe[pip_ind[i]]);
            if (userich)
                pip_RICHbest.push_back(part_RICH_best[pip_ind[i]]);
        }
    }
}

void fill_output_vector_pim(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_pim[i].E() > 0)
        {
            p4_pim_px.push_back(p4_pim[i].Px());
            p4_pim_py.push_back(p4_pim[i].Py());
            p4_pim_pz.push_back(p4_pim[i].Pz());
            p4_pim_E.push_back(p4_pim[i].E());
            pim_det.push_back(pim_detect[i]);
            pim_sec.push_back(pim_sector[i]);
            p4_pim_chi2pid.push_back(pim_chi2pid[i]);
            if(fill_coord_pim == true){
              p4_pim_dcx1.push_back(pim_dcx1[i]); 
              p4_pim_dcy1.push_back(pim_dcy1[i]); 
              p4_pim_dcx2.push_back(pim_dcx2[i]); 
              p4_pim_dcy2.push_back(pim_dcy2[i]); 
              p4_pim_dcx3.push_back(pim_dcx3[i]); 
              p4_pim_dcy3.push_back(pim_dcy3[i]);  
            } 
            if(fill_edge_pim == true){
               p4_pim_dcedge1.push_back(pim_dcedge1[i]); 
               p4_pim_dcedge2.push_back(pim_dcedge2[i]); 
               p4_pim_dcedge3.push_back(pim_dcedge3[i]); 
            }  
            p4_pim_ml.push_back(ml_value(pim_ind[i]));
            p4_pim_beta.push_back(pim_beta[i]);
            p4_pim_LTCC_nphe.push_back(part_CC_LTCC_nphe[pim_ind[i]]);
            if (userich)
                pim_RICHbest.push_back(part_RICH_best[pim_ind[i]]);
        }
    }
}

void fill_output_vector_Kplus(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_Kp[i].E() > 0)
        {
            p4_Kp_px.push_back(p4_Kp[i].Px());
            p4_Kp_py.push_back(p4_Kp[i].Py());
            p4_Kp_pz.push_back(p4_Kp[i].Pz());
            p4_Kp_E.push_back(p4_Kp[i].E());
            Kp_det.push_back(Kp_detect[i]);
            Kp_sec.push_back(Kp_sector[i]);
            p4_Kp_chi2pid.push_back(Kp_chi2pid[i]);
            p4_Kp_ml.push_back(ml_value(Kp_ind[i]));
            p4_Kp_beta.push_back(Kp_beta[i]);
            p4_Kp_LTCC_nphe.push_back(part_CC_LTCC_nphe[Kp_ind[i]]);
            if (userich)
                Kp_RICHbest.push_back(part_RICH_best[Kp_ind[i]]);
        }
    }
}

void fill_output_vector_Kminus(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_Km[i].E() > 0)
        {
            p4_Km_px.push_back(p4_Km[i].Px());
            p4_Km_py.push_back(p4_Km[i].Py());
            p4_Km_pz.push_back(p4_Km[i].Pz());
            p4_Km_E.push_back(p4_Km[i].E());
            Km_det.push_back(Km_detect[i]);
            Km_sec.push_back(Km_sector[i]);
            p4_Km_chi2pid.push_back(Km_chi2pid[i]);
            p4_Km_ml.push_back(ml_value(Km_ind[i]));
            p4_Km_beta.push_back(Km_beta[i]);
            p4_Km_LTCC_nphe.push_back(part_CC_LTCC_nphe[Km_ind[i]]);
            if (userich)
                Km_RICHbest.push_back(part_RICH_best[Km_ind[i]]);
        }
    }
}

void fill_output_vector_photon(void)
{
    for (int i = 0; i < BUFFER; i++)
    {
        if (p4_phot[i].E() > 0)
        {
            p4_phot_px.push_back(p4_phot[i].Px());
            p4_phot_py.push_back(p4_phot[i].Py());
            p4_phot_pz.push_back(p4_phot[i].Pz());
            p4_phot_E.push_back(p4_phot[i].E());
            phot_det.push_back(phot_detect[i]);
            phot_sec.push_back(phot_sector[i]);
            if(fill_coord_phot == true){ 
              p4_phot_pcalv.push_back(phot_pcalv[i]); 
              p4_phot_pcalw.push_back(phot_pcalw[i]);
            } 
        }
    }
}

/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Additional functions:
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double GetTheta(int j)
{
    return acos(part_pz[j] / part_p[j]);
}

double GetPhi(int j)
{
    return 150 * Pival / 180 + atan2(part_py[j] / part_p[j], part_px[j] / part_p[j]);
}

TVector3 GetUVWVector(int j)
{

    double u, v, w, xi, yi, zi;
    double EC_the = 0.4363323;
    double EC_phi;
    double ylow = -182.974;
    double yhi = 189.956;
    double tgrho = 1.95325;
    double sinrho = 0.8901256;
    double cosrho = 0.455715;
    double rot[3][3];

    double x = part_Cal_PCAL_x[j];
    double y = part_Cal_PCAL_y[j];
    double z = part_Cal_PCAL_z[j];

    double at = atan2(y, x);
    if (at < 0)
        at += 2 * Pival;

    double phi = at * 180 / Pival;
    phi = phi + 30.;
    if (phi >= 360.)
        phi = phi - 360.;

    EC_phi = (int)(phi / 60.) * 1.0471975;

    rot[0][0] = cos(EC_the) * cos(EC_phi);
    rot[0][1] = -sin(EC_phi);
    rot[0][2] = sin(EC_the) * cos(EC_phi);
    rot[1][0] = cos(EC_the) * sin(EC_phi);
    rot[1][1] = cos(EC_phi);
    rot[1][2] = sin(EC_the) * sin(EC_phi);
    rot[2][0] = -sin(EC_the);
    rot[2][1] = 0.;
    rot[2][2] = cos(EC_the);

    yi = x * rot[0][0] + y * rot[1][0] + z * rot[2][0];
    xi = x * rot[0][1] + y * rot[1][1] + z * rot[2][1];
    zi = x * rot[0][2] + y * rot[1][2] + z * rot[2][2];

    zi = zi - 510.32;

    u = (yi - ylow) / sinrho;
    v = (yhi - ylow) / tgrho - xi + (yhi - yi) / tgrho;
    w = ((yhi - ylow) / tgrho + xi + (yhi - yi) / tgrho) / 2. / cosrho;

    TVector3 uvw(u, v, w);
    return uvw;
}

double GetTOFmass2(int j, int run)
{
    double tofmass2 = 0;
    if (Beta_charged(j, run) > 0)
        tofmass2 = pow(part_p[j], 2) * (1 - pow(Beta_charged(j, run), 2)) / pow(Beta_charged(j, run), 2);
    return tofmass2;
}

double GetTOFmass2_CD(int j, int run)
{
    double tofmass2 = 0;
    if (Beta_charged_central(j, run) > 0)
        tofmass2 = pow(part_p[j], 2) * (1 - pow(Beta_charged_central(j, run), 2)) / pow(Beta_charged_central(j, run), 2);
    return tofmass2;
}

double Getdvz(int j)
{
    // find index of electron with highest momentum as best vertex approach
    float mom = 0;
    int eind = -1;
    for (int k = 0; k < BUFFER; k++)
    {
        if (FD_eid_all_check[k])
        {
            if (part_p[k] > mom)
            {
                mom = part_p[k];
                eind = k;
            }
        }
    }

    if (eind == -1)
        return -1000; // if no electron exists

    double dvz = 0;
    dvz = part_vz[eind] - part_vz[j];
    return dvz;
}


double Beta_charged(int j, int run)
{
  return part_beta[j];
}

double Beta_charged_central(int j, int run)
{
  return part_beta[j];
}

double Beta_charged_FT(int j, int run)
{
  return part_beta[j];
}


double Beta_neutral(int j, int run)
{
  return part_beta[j];
}

double Beta_neutral_FT(int j, int run)
{
  return part_beta[j];
}



int determineSector(int i)
{
    double phi = 180 / Pival * atan2(part_DC_c2y[i] / sqrt(pow(part_DC_c2x[i], 2) + pow(part_DC_c2y[i], 2) + pow(part_DC_c2z[i], 2)), part_DC_c2x[i] / sqrt(pow(part_DC_c2x[i], 2) + pow(part_DC_c2y[i], 2) + pow(part_DC_c2z[i], 2)));

    if (phi < 30 && phi >= -30)
    {
        return 1;
    }
    else if (phi < 90 && phi >= 30)
    {
        return 2;
    }
    else if (phi < 150 && phi >= 90)
    {
        return 3;
    }
    else if (phi >= 150 || phi < -150)
    {
        return 4;
    }
    else if (phi < -90 && phi >= -150)
    {
        return 5;
    }
    else if (phi < -30 && phi >= -90)
    {
        return 6;
    }

    return 0;
}


bool good_sc_paddle(int j){    // reject bad paddles

  int component_layer1 = part_FTOF_component_layer1[j]; 
  int component_layer2 = part_FTOF_component_layer2[j]; 
  int component_layer3 = part_FTOF_component_layer3[j]; 

  int sector_layer1 = part_FTOF_sector_layer1[j];
  int sector_layer2 = part_FTOF_sector_layer2[j];
  int sector_layer3 = part_FTOF_sector_layer3[j];

  int component_CTOF = part_CTOF_component[j];

  if(sector_layer1 == 2 && component_layer1 == 6)   return false;
  if(sector_layer1 == 2 && component_layer1 == 10)  return false;
  if(sector_layer1 == 2 && component_layer1 == 19)  return false;
  if(sector_layer1 == 6 && component_layer1 == 21)  return false;

  if(sector_layer3 == 1 && component_layer3 == 5)   return false;
  if(sector_layer3 == 2 && component_layer3 == 5)   return false;
  if(sector_layer3 == 6 && component_layer3 == 3)   return false;

  if(component_CTOF == 14)                          return false;

  return true;

}



double ml_value(int i)
{
/*
    float part_beta_buffer, part_chi2pid_buffer;
    float part_p_buffer;
    float part_pid_buffer;

    float part_Cal_PCAL_energy_buffer, part_Cal_ECin_energy_buffer, part_Cal_ECout_energy_buffer, part_Cal_energy_total_buffer;
    float part_Cal_PCAL_time_buffer, part_Cal_ECin_time_buffer, part_Cal_ECout_time_buffer;
    float part_Cal_PCAL_path_buffer, part_Cal_ECin_path_buffer, part_Cal_ECout_path_buffer;

    float part_Cal_PCAL_m2u_buffer, part_Cal_ECin_m2u_buffer, part_Cal_ECout_m2u_buffer;
    float part_Cal_PCAL_m2v_buffer, part_Cal_ECin_m2v_buffer, part_Cal_ECout_m2v_buffer;
    float part_Cal_PCAL_m2w_buffer, part_Cal_ECin_m2w_buffer, part_Cal_ECout_m2w_buffer;
    float part_Cal_PCAL_m3u_buffer, part_Cal_ECin_m3u_buffer, part_Cal_ECout_m3u_buffer;
    float part_Cal_PCAL_m3v_buffer, part_Cal_ECin_m3v_buffer, part_Cal_ECout_m3v_buffer;
    float part_Cal_PCAL_m3w_buffer, part_Cal_ECin_m3w_buffer, part_Cal_ECout_m3w_buffer;
    float part_Cal_PCAL_status_buffer, part_Cal_ECin_status_buffer, part_Cal_ECout_status_buffer;

    int part_CC_HTCC_sector_buffer, part_CC_HTCC_nphe_buffer;
    float part_CC_HTCC_time_buffer, part_CC_HTCC_path_buffer;
    float part_CC_HTCC_theta_buffer, part_CC_HTCC_phi_buffer;

    int part_CC_LTCC_sector_buffer, part_CC_LTCC_nphe_buffer;
    float part_CC_LTCC_time_buffer, part_CC_LTCC_path_buffer;
    float part_CC_LTCC_theta_buffer, part_CC_LTCC_phi_buffer;

    int part_FTOF_layer_buffer;
    int part_FTOF_sector_layer1_buffer, part_FTOF_sector_layer2_buffer, part_FTOF_sector_layer3_buffer;
    int part_FTOF_component_layer1_buffer, part_FTOF_component_layer2_buffer, part_FTOF_component_layer3_buffer;
    float part_FTOF_energy_buffer, part_FTOF_time_buffer, part_FTOF_path_buffer;
    float part_FTOF_energy_layer1_buffer, part_FTOF_time_layer1_buffer, part_FTOF_path_layer1_buffer;
    float part_FTOF_energy_layer3_buffer, part_FTOF_time_layer3_buffer, part_FTOF_path_layer3_buffer;

    float part_theta_buffer;

    float part_Cal_PCAL_velocity_buffer;

    float part_Cal_ECin_velocity_buffer;

    float part_Cal_ECout_velocity_buffer;

    float part_CC_HTCC_nphe_float_buffer;

    float part_CC_HTCC_velocity_buffer;

    float part_chi2trk_buffer;

    float part_FTOF_velocity_buffer;

    float part_FTOF_velocity_layer1_buffer;

    float part_FTOF_velocity_layer3_buffer;

    float partLUND_pid_buffer;

    part_p_buffer = part_p[i];
    part_beta_buffer = part_beta[i];
    part_pid_buffer = part_pid[i];
    part_chi2pid_buffer = part_chi2pid[i];

    part_Cal_PCAL_energy_buffer = part_Cal_PCAL_energy[i];
    part_Cal_ECin_energy_buffer = part_Cal_ECin_energy[i];
    part_Cal_ECout_energy_buffer = part_Cal_ECout_energy[i];
    part_Cal_energy_total_buffer = part_Cal_energy_total[i];
    part_Cal_PCAL_time_buffer = part_Cal_PCAL_time[i];
    part_Cal_ECin_time_buffer = part_Cal_ECin_time[i];
    part_Cal_ECout_time_buffer = part_Cal_ECout_time[i];
    part_Cal_PCAL_path_buffer = part_Cal_PCAL_path[i];
    part_Cal_ECin_path_buffer = part_Cal_ECin_path[i];
    part_Cal_ECout_path_buffer = part_Cal_ECout_path[i];

    part_Cal_PCAL_m2u_buffer = part_Cal_PCAL_m2u[i];
    part_Cal_PCAL_m2v_buffer = part_Cal_PCAL_m2v[i];
    part_Cal_PCAL_m2w_buffer = part_Cal_PCAL_m2w[i];
    part_Cal_PCAL_m3u_buffer = part_Cal_PCAL_m3u[i];
    part_Cal_PCAL_m3v_buffer = part_Cal_PCAL_m3v[i];
    part_Cal_PCAL_m3w_buffer = part_Cal_PCAL_m3w[i];
    part_Cal_ECin_m2u_buffer = part_Cal_ECin_m2u[i];
    part_Cal_ECin_m2v_buffer = part_Cal_ECin_m2v[i];
    part_Cal_ECin_m2w_buffer = part_Cal_ECin_m2w[i];
    part_Cal_ECin_m3u_buffer = part_Cal_ECin_m3u[i];
    part_Cal_ECin_m3v_buffer = part_Cal_ECin_m3v[i];
    part_Cal_ECin_m3w_buffer = part_Cal_ECin_m3w[i];
    part_Cal_ECout_m2u_buffer = part_Cal_ECout_m2u[i];
    part_Cal_ECout_m2v_buffer = part_Cal_ECout_m2v[i];
    part_Cal_ECout_m2w_buffer = part_Cal_ECout_m2w[i];
    part_Cal_ECout_m3u_buffer = part_Cal_ECout_m3u[i];
    part_Cal_ECout_m3v_buffer = part_Cal_ECout_m3v[i];
    part_Cal_ECout_m3w_buffer = part_Cal_ECout_m3w[i];

    part_CC_HTCC_nphe_buffer = part_CC_HTCC_nphe[i];
    part_CC_HTCC_time_buffer = part_CC_HTCC_time[i];
    part_CC_HTCC_path_buffer = part_CC_HTCC_path[i];

    part_FTOF_energy_buffer = part_FTOF_energy[i];
    part_FTOF_time_buffer = part_FTOF_time[i];
    part_FTOF_path_buffer = part_FTOF_path[i];
    part_FTOF_energy_layer1_buffer = part_FTOF_energy_layer1[i];
    part_FTOF_time_layer1_buffer = part_FTOF_time_layer1[i];
    part_FTOF_path_layer1_buffer = part_FTOF_path_layer1[i];
    part_FTOF_energy_layer3_buffer = part_FTOF_energy_layer3[i];
    part_FTOF_time_layer3_buffer = part_FTOF_time_layer3[i];
    part_FTOF_path_layer3_buffer = part_FTOF_path_layer3[i];
    part_FTOF_layer_buffer = part_FTOF_layer[i];

    part_theta_buffer = part_theta[i];

    part_Cal_PCAL_velocity_buffer = part_Cal_PCAL_path_buffer / part_Cal_PCAL_time_buffer;

    if (std::isnan(part_Cal_PCAL_velocity_buffer))
        part_Cal_PCAL_velocity_buffer = 0.;

    part_Cal_ECin_velocity_buffer = part_Cal_ECin_path_buffer / part_Cal_ECin_time_buffer;

    if (std::isnan(part_Cal_ECin_velocity_buffer))
        part_Cal_ECin_velocity_buffer = 0.;

    part_Cal_ECout_velocity_buffer = part_Cal_ECout_path_buffer / part_Cal_ECout_time_buffer;

    if (std::isnan(part_Cal_ECout_velocity_buffer))
        part_Cal_ECout_velocity_buffer = 0.;

    part_CC_HTCC_nphe_float_buffer = (float)part_CC_HTCC_nphe_buffer;

    part_CC_HTCC_velocity_buffer = part_CC_HTCC_path_buffer / part_CC_HTCC_time_buffer;

    part_FTOF_velocity_buffer = part_FTOF_path_buffer / part_FTOF_time_buffer;

    part_FTOF_velocity_layer1_buffer = part_FTOF_path_layer1_buffer / part_FTOF_time_layer1_buffer;

    part_FTOF_velocity_layer3_buffer = part_FTOF_path_layer3_buffer / part_FTOF_time_layer3_buffer;

    if (std::isnan(part_CC_HTCC_velocity_buffer))
        part_CC_HTCC_velocity_buffer = 0.;

    if (std::isnan(part_FTOF_velocity_buffer))
        part_FTOF_velocity_buffer = 0.;

    if (std::isnan(part_FTOF_velocity_layer1_buffer))
        part_FTOF_velocity_layer1_buffer = 0.;

    if (std::isnan(part_FTOF_velocity_layer3_buffer))
        part_FTOF_velocity_layer3_buffer = 0.;

    partLUND_pid_buffer = 0.0;

    TMVA::Reader *reader = new TMVA::Reader("!Color:Silent");

    reader->AddVariable("part_beta", &part_beta_buffer);
    reader->AddVariable("part_p", &part_p_buffer);

    reader->AddVariable("part_Cal_PCAL_energy", &part_Cal_PCAL_energy_buffer);
    reader->AddVariable("part_Cal_ECin_energy", &part_Cal_ECin_energy_buffer);
    reader->AddVariable("part_Cal_ECout_energy", &part_Cal_ECout_energy_buffer);
    reader->AddVariable("part_Cal_energy_total", &part_Cal_energy_total_buffer);
    reader->AddVariable("part_Cal_PCAL_path/part_Cal_PCAL_time", &part_Cal_PCAL_velocity_buffer);
    reader->AddVariable("part_Cal_ECin_path/part_Cal_ECin_time", &part_Cal_ECin_velocity_buffer);
    reader->AddVariable("part_Cal_ECout_path/part_Cal_ECout_time", &part_Cal_ECout_velocity_buffer);

    reader->AddVariable("part_Cal_PCAL_m2u", &part_Cal_PCAL_m2u_buffer);
    reader->AddVariable("part_Cal_ECin_m2u", &part_Cal_ECin_m2u_buffer);
    reader->AddVariable("part_Cal_ECout_m2u", &part_Cal_ECout_m2u_buffer);
    reader->AddVariable("part_Cal_PCAL_m2v", &part_Cal_PCAL_m2v_buffer);
    reader->AddVariable("part_Cal_ECin_m2v", &part_Cal_ECin_m2v_buffer);
    reader->AddVariable("part_Cal_ECout_m2v", &part_Cal_ECout_m2v_buffer);
    reader->AddVariable("part_Cal_PCAL_m2w", &part_Cal_PCAL_m2w_buffer);
    reader->AddVariable("part_Cal_ECin_m2w", &part_Cal_ECin_m2w_buffer);
    reader->AddVariable("part_Cal_ECout_m2w", &part_Cal_ECout_m2w_buffer);
    reader->AddVariable("part_Cal_PCAL_m3u", &part_Cal_PCAL_m3u_buffer);
    reader->AddVariable("part_Cal_ECin_m3u", &part_Cal_ECin_m3u_buffer);
    reader->AddVariable("part_Cal_ECout_m3u", &part_Cal_ECout_m3u_buffer);
    reader->AddVariable("part_Cal_PCAL_m3v", &part_Cal_PCAL_m3v_buffer);
    reader->AddVariable("part_Cal_ECin_m3v", &part_Cal_ECin_m3v_buffer);
    reader->AddVariable("part_Cal_ECout_m3v", &part_Cal_ECout_m3v_buffer);
    reader->AddVariable("part_Cal_PCAL_m3w", &part_Cal_PCAL_m3w_buffer);
    reader->AddVariable("part_Cal_ECin_m3w", &part_Cal_ECin_m3w_buffer);
    reader->AddVariable("part_Cal_ECout_m3w", &part_Cal_ECout_m3w_buffer);

    reader->AddVariable("part_CC_HTCC_nphe", &part_CC_HTCC_nphe_float_buffer);
    reader->AddVariable("part_CC_HTCC_path/part_CC_HTCC_time", &part_CC_HTCC_velocity_buffer);

    reader->AddVariable("part_chi2pid", &part_chi2pid_buffer);
    reader->AddVariable("part_pid", &part_pid_buffer);

    reader->AddVariable("part_FTOF_energy", &part_FTOF_energy_buffer);
    reader->AddVariable("part_FTOF_path/part_FTOF_time", &part_FTOF_velocity_buffer);
    reader->AddVariable("part_FTOF_energy_layer1", &part_FTOF_energy_layer1_buffer);
    reader->AddVariable("part_FTOF_path_layer1/part_FTOF_time_layer1", &part_FTOF_velocity_layer1_buffer);
    reader->AddVariable("part_FTOF_energy_layer3", &part_FTOF_energy_layer3_buffer);
    reader->AddVariable("part_FTOF_path_layer3/part_FTOF_time_layer3", &part_FTOF_velocity_layer3_buffer);

    reader->AddSpectator("partLUND_pid", &partLUND_pid_buffer);

    if (part_charge[i] > 0)
    {
        reader->BookMVA("DNN", "./dataset/weight3FinValFTOF/TMVAClassification_DNN.weights.xml");
    }
    else
    {
        reader->BookMVA("DNN", "/home/aron/CLAS/analysis/dataset/home/aron/CLAS/analysis/dataset/Km/TMVAClassification_DNN.weights.xml");
    }

    if (!std::isnan(part_Cal_PCAL_energy_buffer) && !std::isnan(part_Cal_ECin_energy_buffer) && !std::isnan(part_Cal_ECout_energy_buffer) && !std::isnan(part_Cal_energy_total_buffer) && !std::isnan(part_Cal_PCAL_time_buffer) && !std::isnan(part_Cal_ECin_time_buffer) && !std::isnan(part_Cal_ECout_time_buffer) && !std::isnan(part_Cal_PCAL_path_buffer) && !std::isnan(part_Cal_ECin_path_buffer) && !std::isnan(part_Cal_ECout_path_buffer) && !std::isnan(part_Cal_PCAL_m2u_buffer) && !std::isnan(part_Cal_ECin_m2u_buffer) && !std::isnan(part_Cal_ECout_m2u_buffer) && !std::isnan(part_Cal_PCAL_m2v_buffer) && !std::isnan(part_Cal_ECin_m2v_buffer) && !std::isnan(part_Cal_ECout_m2v_buffer) && !std::isnan(part_Cal_PCAL_m2w_buffer) && !std::isnan(part_Cal_ECin_m2w_buffer) && !std::isnan(part_Cal_ECout_m2w_buffer) && !std::isnan(part_Cal_PCAL_m3u_buffer) && !std::isnan(part_Cal_ECin_m3u_buffer) && !std::isnan(part_Cal_ECout_m3u_buffer) && !std::isnan(part_Cal_PCAL_m3v_buffer) && !std::isnan(part_Cal_ECin_m3v_buffer) && !std::isnan(part_Cal_ECout_m3v_buffer) && !std::isnan(part_Cal_PCAL_m3w_buffer) && !std::isnan(part_Cal_ECin_m3w_buffer) && !std::isnan(part_Cal_ECout_m3w_buffer) && !std::isnan(part_CC_HTCC_nphe_buffer) && !std::isnan(part_CC_HTCC_time_buffer) && !std::isnan(part_CC_HTCC_path_buffer) && !std::isnan(part_beta_buffer) && !std::isnan(part_p_buffer) && part_theta_buffer > 0.09 && part_theta_buffer < 0.61 && part_beta_buffer < 1.1 && part_beta_buffer > 0. && part_Cal_PCAL_energy_buffer >= 0. && part_Cal_ECin_energy_buffer >= 0. && part_Cal_ECout_energy_buffer >= 0. && part_Cal_energy_total_buffer < 0.4 && part_Cal_PCAL_time_buffer >= 0. && part_Cal_ECin_time_buffer >= 0. && part_Cal_ECout_time_buffer >= 0. && part_Cal_PCAL_path_buffer >= 0. && part_Cal_ECin_path_buffer >= 0. && part_Cal_ECout_path_buffer >= 0. && part_Cal_PCAL_m2u_buffer >= 0. && part_Cal_ECin_m2u_buffer >= 0. && part_Cal_ECout_m2u_buffer >= 0. && part_Cal_PCAL_m2w_buffer >= 0. && part_Cal_ECin_m2w_buffer >= 0. && part_Cal_ECout_m2w_buffer >= 0. && part_Cal_PCAL_m2v_buffer >= 0. && part_Cal_ECin_m2v_buffer >= 0. && part_Cal_ECout_m2v_buffer >= 0. && part_CC_HTCC_nphe_buffer >= 0. && part_CC_HTCC_time_buffer >= 0. && part_CC_HTCC_path_buffer >= 0.)
    {
        const double return_value = reader->EvaluateMVA("DNN");
        delete reader;
        return return_value;
    }

    delete reader;
*/
    return -99;
}


/*
double part_Cal_PCAL_velocity = part_Cal_PCAL_path[i] / part_Cal_PCAL_time[i];

if (std::isnan(part_Cal_PCAL_velocity))
    part_Cal_PCAL_velocity = 0.;

double part_Cal_ECin_velocity = part_Cal_ECin_path[i] / part_Cal_ECin_time[i];

if (std::isnan(part_Cal_ECin_velocity))
    part_Cal_ECin_velocity = 0.;

double part_Cal_ECout_velocity = part_Cal_ECout_path[i] / part_Cal_ECout_time[i];

if (std::isnan(part_Cal_ECout_velocity))
    part_Cal_ECout_velocity = 0.;

double part_CC_HTCC_nphe_float = static_cast<double>(part_CC_HTCC_nphe[i]);

double part_CC_HTCC_velocity = part_CC_HTCC_path[i] / part_CC_HTCC_time[i];

double part_FTOF_velocity = part_FTOF_path[i] / part_FTOF_time[i];

double part_FTOF_velocity_layer1 = part_FTOF_path_layer1[i] / part_FTOF_time_layer1[i];

double part_FTOF_velocity_layer3 = part_FTOF_path_layer3[i] / part_FTOF_time_layer3[i];

double part_pid_float = static_cast<double>(part_pid[i]);

if (std::isnan(part_CC_HTCC_velocity))
    part_CC_HTCC_velocity = 0.;

if (std::isnan(part_FTOF_velocity))
    part_FTOF_velocity = 0.;

if (std::isnan(part_FTOF_velocity_layer1))
    part_FTOF_velocity_layer1 = 0.;

if (std::isnan(part_FTOF_velocity_layer3))
    part_FTOF_velocity_layer3 = 0.;

if (!std::isnan(part_Cal_PCAL_energy[i]) && !std::isnan(part_Cal_ECin_energy[i]) && !std::isnan(part_Cal_ECout_energy[i]) && !std::isnan(part_Cal_energy_total[i]) && !std::isnan(part_Cal_PCAL_time[i]) && !std::isnan(part_Cal_ECin_time[i]) && !std::isnan(part_Cal_ECout_time[i]) && !std::isnan(part_Cal_PCAL_path[i]) && !std::isnan(part_Cal_ECin_path[i]) && !std::isnan(part_Cal_ECout_path[i]) && !std::isnan(part_Cal_PCAL_m2u[i]) && !std::isnan(part_Cal_ECin_m2u[i]) && !std::isnan(part_Cal_ECout_m2u[i]) && !std::isnan(part_Cal_PCAL_m2v[i]) && !std::isnan(part_Cal_ECin_m2v[i]) && !std::isnan(part_Cal_ECout_m2v[i]) && !std::isnan(part_Cal_PCAL_m2w[i]) && !std::isnan(part_Cal_ECin_m2w[i]) && !std::isnan(part_Cal_ECout_m2w[i]) && !std::isnan(part_Cal_PCAL_m3u[i]) && !std::isnan(part_Cal_ECin_m3u[i]) && !std::isnan(part_Cal_ECout_m3u[i]) && !std::isnan(part_Cal_PCAL_m3v[i]) && !std::isnan(part_Cal_ECin_m3v[i]) && !std::isnan(part_Cal_ECout_m3v[i]) && !std::isnan(part_Cal_PCAL_m3w[i]) && !std::isnan(part_Cal_ECin_m3w[i]) && !std::isnan(part_Cal_ECout_m3w[i]) && !std::isnan(part_CC_HTCC_nphe[i]) && !std::isnan(part_CC_HTCC_time[i]) && !std::isnan(part_CC_HTCC_path[i]) && !std::isnan(part_beta[i]) && !std::isnan(part_p[i]) && part_theta[i] > 0.09 && part_theta[i] < 0.61 && part_beta[i] < 1.1 && part_beta[i] > 0. && part_Cal_PCAL_energy[i] >= 0. && part_Cal_ECin_energy[i] >= 0. && part_Cal_ECout_energy[i] >= 0. && part_Cal_energy_total[i] < 0.4 && part_Cal_PCAL_time[i] >= 0. && part_Cal_ECin_time[i] >= 0. && part_Cal_ECout_time[i] >= 0. && part_Cal_PCAL_path[i] >= 0. && part_Cal_ECin_path[i] >= 0. && part_Cal_ECout_path[i] >= 0. && part_Cal_PCAL_m2u[i] >= 0. && part_Cal_ECin_m2u[i] >= 0. && part_Cal_ECout_m2u[i] >= 0. && part_Cal_PCAL_m2w[i] >= 0. && part_Cal_ECin_m2w[i] >= 0. && part_Cal_ECout_m2w[i] >= 0. && part_Cal_PCAL_m2v[i] >= 0. && part_Cal_ECin_m2v[i] >= 0. && part_Cal_ECout_m2v[i] >= 0. && part_CC_HTCC_nphe[i] >= 0. && part_CC_HTCC_time[i] >= 0. && part_CC_HTCC_path[i] >= 0.)
{
    std::vector<std::string> inputVars{"part_beta", "part_p", "part_Cal_PCAL_energy", "part_Cal_ECin_energy", "part_Cal_ECout_energy", "part_Cal_energy_total", "part_Cal_PCAL_path/part_Cal_PCAL_time", "part_Cal_ECin_path/part_Cal_ECin_time", "part_Cal_ECout_path/part_Cal_ECout_time", "part_Cal_PCAL_m2u", "part_Cal_ECin_m2u", "part_Cal_ECout_m2u", "part_Cal_PCAL_m2v", "part_Cal_ECin_m2v", "part_Cal_ECout_m2v", "part_Cal_PCAL_m2w", "part_Cal_ECin_m2w", "part_Cal_ECout_m2w", "part_Cal_PCAL_m3u", "part_Cal_ECin_m3u", "part_Cal_ECout_m3u", "part_Cal_PCAL_m3v", "part_Cal_ECin_m3v", "part_Cal_ECout_m3v", "part_Cal_PCAL_m3w", "part_Cal_ECin_m3w", "part_Cal_ECout_m3w", "part_CC_HTCC_nphe", "part_CC_HTCC_path/part_CC_HTCC_time", "part_chi2pid", "part_pid", "part_FTOF_energy", "part_FTOF_path/part_FTOF_time", "part_FTOF_energy_layer1", "part_FTOF_path_layer1/part_FTOF_time_layer1", "part_FTOF_energy_layer3", "part_FTOF_path_layer3/part_FTOF_time_layer3"};

    const ReadDNN *DNNResponse = new ReadDNN(inputVars);

    const std::vector<double> inputVec{part_beta[i], part_p[i], part_Cal_PCAL_energy[i], part_Cal_ECin_energy[i], part_Cal_ECout_energy[i], part_Cal_energy_total[i], part_Cal_PCAL_velocity, part_Cal_ECin_velocity, part_Cal_ECout_velocity, part_Cal_PCAL_m2u[i], part_Cal_ECin_m2u[i], part_Cal_ECout_m2u[i], part_Cal_PCAL_m2v[i], part_Cal_ECin_m2v[i], part_Cal_ECout_m2v[i], part_Cal_PCAL_m2w[i], part_Cal_ECin_m2w[i], part_Cal_ECout_m2w[i], part_Cal_PCAL_m3u[i], part_Cal_ECin_m3u[i], part_Cal_ECout_m3u[i], part_Cal_PCAL_m3v[i], part_Cal_ECin_m3v[i], part_Cal_ECout_m3v[i], part_Cal_PCAL_m3w[i], part_Cal_ECin_m3w[i], part_Cal_ECout_m3w[i], part_CC_HTCC_nphe_float, part_CC_HTCC_velocity, part_chi2pid[i], part_pid_float, part_FTOF_energy[i], part_FTOF_velocity, part_FTOF_energy_layer1[i], part_FTOF_velocity_layer1, part_FTOF_energy_layer3[i], part_FTOF_velocity_layer3};

    const double return_value = DNNResponse->GetMvaValue(inputVec);

    delete DNNResponse;

    return return_value;
}

return -99;


*/
