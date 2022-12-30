#ifndef CONFIGURATION__H
#define CONFIGURATION__H

#include "util/globalMacros.h"
#include "util/globalTypesClasses.h"
#include "util/OneDimensionalElement.h"

// weights of the star values for one side of the interface
// f_ -> factor of
// w: stands for v or u depending on the formulation
class StarW1s {
public:
  // sigma_nStar = ss_f_sigmaL * sigmaL + ss_f_sigmaR * sigmaR + ss_f_wL * wL +
  // ss_f_wR * wR
  NUMBR ss_f_sigmaL, ss_f_sigmaR, ss_f_wL, ss_f_wR;
  // w_Star = ws_f_sigmaL * sigmaL + ws_f_sigmaR * sigmaR + ws_f_wL * wL +
  // ws_f_wR * wR
  NUMBR ws_f_sigmaL, ws_f_sigmaR, ws_f_wL, ws_f_wR;
};

class StarW2s {
public:
  // side_weights[SDL]: the weights of the star values for the left side
  // side_weights[SDR]: // right side
  StarW1s side_weights[NUM_SIDES];
};

class Configuration {
public:
  Configuration();
  // from the main inputs -> compute data in "computed" block
  void Initialize_Configuration();

  // inputs
  unsigned int polyOrder;
  // whether we want the mass-lumpting option or not
  bool lumpMass;

  WeakFormulationType wf_type;
  // star options for DG method
  StarOption sOption;
  DG_epsilonT dg_eps;
  // for DG, non-Riemann star options sigma*n has a penalty multiplying (wR -
  // wL), that weight is computed from eta below etaI: for interior interfaces
  // etaB: for boundary faces
  double etaI, etaB;

  // computed from inputs
  bool isDG;
  unsigned int num_fields;
  // position of U in terms of unknowns of an element -> for U it's always 0
  unsigned int field_Start_U;
  // field_Start_V: for 1F formulation V is computed from U -> position in
  // element unknowns is 0 for 2F -> V is position 1 within the fields of the
  // elements
  unsigned int field_Start_V;

  // w* = v* -> isW_velocity = true, w* = u* -> isW_velocity = false
  bool isW_velocity;

  // vHat (rho vDot + damping v - div. S - source) -> weight_BLM_is_velocity =
  // true uHat // -> weight_BLM_is_velocity = false
  bool weight_BLM_is_velocity;
  // what is the position of the weight relative to unknown fields of the
  // element
  unsigned int field_pos_weight_BLM;
  bool isBlochModeAnalysis;
  Dcomplex gamma, gamma_inv; // gamma = exp(i k L) // k wavenumber, L is domain
                             // length -> needed for Bloch analysis

  bool isHelmholtz;
  double omega;

  /// large class members
  // we read / have already initialized [he, E, rho, damping] of elements[i].elementProps for each element and calculate Z, c, etc. from this
  vector<OneDimensionalElement> elements;
  unsigned int num_elements;

  ////////////////////////////////////////////
  // Computed
  OneDimensionalParentElement parentElement;
  unsigned int ndof_parent_element;
  unsigned int ndof_element;
  unsigned int ndof_domain;
  // I think we can remove this for 1D DG/CFEM
  vector<vector< int > > edof_2_globalDofMap;
	// these are indexed by ei (element index) and store the start and end dof of element (element goes grom  = start; < end)
  vector<int> element_start_dof;
  vector<int> element_end_dof;


  /////////////////////////////////// Global, K, M, C
  // stores whether we need a damping matrix
  bool b_hasNonZeroDampingMatrix;

  ///	Auxiliary functions

  // inputs: left_ep, right_ep -> const makes it clearer that these are inputs
  //		   bool insideDomainInterface: periodic and Bloch BCs are
  // interior
  // interfaces. Bloch BC involves the value gamma 				other
  // interfaces that fall inside the domain are labeled as
  // "insideDomainInterface = true" for Bloch and periodic BC
  // insideDomainInterface = false outputs: twoSideWeights
  void Compute_DG_Star_Weights_4_Inteior_Interface(
      const ElementProperties &left_ep, const ElementProperties &right_ep,
      bool insideDomainInterface,
      StarW2s &twoSideWeights) const; // const here says that this object
                                      // (configuration) is not changing itself

private:
	// for 1F-formulation we end up with Ma'' + Ca' + Ka = F -> damping only if we have nonzero damping at element level or we have 1Dvs formulation
	// for 2F-formulation Ma' + Ka = F -> no damping
	bool HasNonZeroDampingMatrix() const;
	void Initialize_Elements();
};

#endif