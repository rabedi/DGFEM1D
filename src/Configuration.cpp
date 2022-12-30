#include "Configuration.h"

Configuration::Configuration() {
  wf_type = DG_2FUV;
  sOption = so_Riemann;
  dg_eps = dg_eps_p1;
  etaI = 1.0;
  etaB = 1.0;
  Initialize_Configuration();
}

void Configuration::Initialize_Configuration() {
  isDG = (wf_type != cfem_1D);
  weight_BLM_is_velocity = (wf_type == DG_2FUV);
  num_fields = 1;
  field_Start_U = 0;
  field_Start_V = 0;
  field_pos_weight_BLM = 0;
  if (wf_type == DG_2FUV) {
    num_fields = 2;
    field_pos_weight_BLM = 1;
  }
}

void Configuration::Compute_DG_Star_Weights_4_Inteior_Interface(
    const ElementProperties &left_ep, const ElementProperties &right_ep,
    bool insideDomainInterface, StarW2s &twoSideWeights) const {
  // calculates sigmaStar (not sigmaStar_n = sigmaStar . n) and wStar weight
  StarW1s shared_sigmaStar_wStarWeights;
  if (sOption == so_Riemann) {
    double Zl = left_ep.Z, Zr = right_ep.Z;
    double ZlpZr_inv = 1.0 / (Zl + Zr);
    double Zl_div_ZlpZr = Zl * Zl_div_ZlpZr, Zr_div_ZlpZr = Zr * Zl_div_ZlpZr;

    // sigma* weights
    shared_sigmaStar_wStarWeights.ss_f_sigmaL = Zl_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ss_f_sigmaR = Zr_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ss_f_wR = Zl * Zr * ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ss_f_wL =
        -shared_sigmaStar_wStarWeights.ss_f_wR;

    // Velocity* weights (still not necessarily w weights)
    shared_sigmaStar_wStarWeights.ws_f_sigmaL = -ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR = ZlpZr_inv;
    shared_sigmaStar_wStarWeights.ws_f_wR = Zl_div_ZlpZr;
    shared_sigmaStar_wStarWeights.ws_f_wL = Zr_div_ZlpZr;
  } else {
    // taking care of the penality term
    if (etaI > 1e-15) {
      double eta = 0.0; // penality term - absolute value
      if (isW_velocity) {
        double Zl = left_ep.Z, Zr = right_ep.Z;
        double ZlpZr_inv = 1.0 / (Zl + Zr);
        eta = etaI * Zl * Zr * ZlpZr_inv;
      } else // w is u -> elliptic start
      {
        double ETilde = 0.5 * (left_ep.E + right_ep.E);
        double hTilde = 0.5 * (left_ep.hE + right_ep.hE);
        eta = etaI * ETilde / hTilde;
      }
      // sigma* = eta * (wR - wL)
      shared_sigmaStar_wStarWeights.ss_f_wR = -eta;
      shared_sigmaStar_wStarWeights.ss_f_wL = eta;
    }
    shared_sigmaStar_wStarWeights.ws_f_sigmaL = 0.0;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR = 0.0;

    if (sOption == so_Central) {
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.5;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.5;

      shared_sigmaStar_wStarWeights.ws_f_wR = 0.5;
      shared_sigmaStar_wStarWeights.ws_f_wL = 0.5;
    } else if (sOption == so_Alternating_sL) {
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 1.0;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 0.0;

      shared_sigmaStar_wStarWeights.ws_f_wR = 0.0;
      shared_sigmaStar_wStarWeights.ws_f_wL = 1.0;
    } 
	else if (sOption == so_Alternating_sR) 
	{
      shared_sigmaStar_wStarWeights.ss_f_sigmaL = 0.0;
      shared_sigmaStar_wStarWeights.ss_f_sigmaR = 1.0;

      shared_sigmaStar_wStarWeights.ws_f_wR = 1.0;
      shared_sigmaStar_wStarWeights.ws_f_wL = 0.0;
    } else {
      cout << "sOption\t" << sOption << '\n';
      THROW("Invalid sOption\n");
    }
  }

  /// if w is u (Helmholtz), we need to use i omega to adjust the weights
  if (isHelmholtz) {
    Dcomplex iomega = omega * Icomp, iomega_inv = 1.0 / iomega;

    shared_sigmaStar_wStarWeights.ss_f_wR *= iomega;
    shared_sigmaStar_wStarWeights.ss_f_wL *= iomega;

    // Velocity* weights (still not necessarily w weights)
    shared_sigmaStar_wStarWeights.ws_f_sigmaL *= iomega_inv;
    shared_sigmaStar_wStarWeights.ws_f_sigmaR *= iomega_inv;
  }
  twoSideWeights.side_weights[SDL] = shared_sigmaStar_wStarWeights;
  twoSideWeights.side_weights[SDR] = shared_sigmaStar_wStarWeights;
  // need to multiply sigma*n parts of the right by -1, since for the right side
  // n = -1
  twoSideWeights.side_weights[SDR].ss_f_sigmaL =
      -twoSideWeights.side_weights[SDR].ss_f_sigmaL;
  twoSideWeights.side_weights[SDR].ss_f_sigmaR =
      -twoSideWeights.side_weights[SDR].ss_f_sigmaR;
  twoSideWeights.side_weights[SDR].ss_f_wL =
      -twoSideWeights.side_weights[SDR].ss_f_wL;
  twoSideWeights.side_weights[SDR].ss_f_wR =
      -twoSideWeights.side_weights[SDR].ss_f_wR;

  if ((!isBlochModeAnalysis) || insideDomainInterface) // ready to return
    return;

  // deal with Bloch mode analysis where gamma is involved
  twoSideWeights.side_weights[SDL].ss_f_sigmaR *= gamma;
  twoSideWeights.side_weights[SDL].ss_f_wR *= gamma;
  twoSideWeights.side_weights[SDL].ws_f_sigmaR *= gamma;
  twoSideWeights.side_weights[SDL].ws_f_wR *= gamma;

  twoSideWeights.side_weights[SDR].ss_f_sigmaL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ss_f_wL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ws_f_sigmaL *= gamma_inv;
  twoSideWeights.side_weights[SDR].ws_f_wL *= gamma_inv;
}

bool Configuration::HasNonZeroDampingMatrix() const
{
	if (wf_type == DG_2FUV)
		return false;
	if (wf_type == DG_1F_vStar)
		return true;

	static double tol_4_damping = 1000 * DBL_MIN;
	if ((wf_type == cfem_1D) || (wf_type == DG_1F_uStar))
	{
		// depends on the material damping
		for (unsigned int ei = 0; ei < num_elements; ++ei)
		{
			if (elements[ei].elementProps.damping > tol_4_damping)
				return true;
		}
		return false;
	}
}

void Configuration::Initialize_Elements()
{
	// initialize the parent
	parentElement.Initialize(polyOrder, lumpMass);
	ndof_parent_element = parentElement.ndof;
	ndof_element = ndof_parent_element * num_fields;

	edof_2_globalDofMap.resize(num_elements);
	element_start_dof.resize(num_elements);
	element_end_dof.resize(num_elements);

	edof_2_globalDofMap.clear();

	if (wf_type != cfem_1D)
	{
		ndof_domain = ndof_element * num_elements;
		vector<int> dof_map_element(ndof_element);
		unsigned cntr = 0;
		for (unsigned int ei = 0; ei < num_elements; ++ei)
		{
			element_start_dof[ei] = cntr;
			element_end_dof[ei] = cntr + ndof_element;
			for (unsigned int j = 0; j < ndof_element; ++j)
				dof_map_element[j] = cntr++;
			edof_2_globalDofMap.push_back(dof_map_element);
		}
	}
	else
	{
		THROW("Add CFEM block later\n");
	}

	b_hasNonZeroDampingMatrix = HasNonZeroDampingMatrix();

	for (unsigned int ei = 0; ei < num_elements; ++ei)
	{
		OneDimensionalElement* ePtr = &elements[ei];
		ePtr->elementProps.Initialize_ElementProperties();

		// these building blocks are mpe, kpe with 1F size but with real element geometry (size) and no material property
		MATRIX meBuildingBlock, keBuildingBlock;
		double Jacobian = 0.5 * ePtr->elementProps.hE;
		keBuildingBlock.Multiply(parentElement.kpe, Jacobian);
		meBuildingBlock.Multiply(parentElement.kpe, Jacobian);


		ePtr->ke.resize(ndof_element);
		ePtr->ke = 0.0;
		ePtr->me.resize(ndof_element);
		ePtr->me = 0.0;
		if (b_hasNonZeroDampingMatrix)
		{
			ePtr->ce.resize(ndof_element);
			ePtr->ce = 0.0;
		}

		if (wf_type == DG_2FUV)
		{
			double alpha = 1.0 / ePtr->elementProps.time_e;
			alpha *= alpha;
			for (unsigned int i = 0; i < parentElement.ndof; i)
			{
				for (unsigned int j = 0; j < parentElement.ndof; ++j)
				{
					/// stiffness term
					// -alpha uHat * v
					ePtr->ke[i][j + parentElement.ndof] = -alpha * meBuildingBlock[i][j];
					// vHat * damping * v
					ePtr->ke[i + parentElement.ndof][j + parentElement.ndof] = ePtr->elementProps.damping * meBuildingBlock[i][j];
					// grad vHat * sigma
					ePtr->ke[i + parentElement.ndof][j] = ePtr->elementProps.E * keBuildingBlock[i][j];

					/// mass term				
					// alpha uHat * uDot
					ePtr->me[i][j] = alpha * meBuildingBlock[i][j];
					// vHat * vDot
					ePtr->me[i + parentElement.ndof][j + parentElement.ndof] = ePtr->elementProps.rho * meBuildingBlock[i][j];
				}
			}
		}
		else //if ((wf_type == DG_1F_vStar) || (wf_type == DG_1F_uStar) || (wf_type == cfem_1D))
		{
			for (unsigned int i = 0; i < parentElement.ndof; i)
			{
				for (unsigned int j = 0; j < parentElement.ndof; ++j)
				{
					/// stiffness term
					// grad uHat * sigma
					ePtr->ke[i][j] = ePtr->elementProps.E * keBuildingBlock[i][j];
					/// mass term				
					// uHat * vDot = uHat * uDDot
					ePtr->me[i][j] = ePtr->elementProps.rho * meBuildingBlock[i][j];

					/// damping term				
					// uHat * damping * v = uHat * damping * uDot
					if (b_hasNonZeroDampingMatrix)
						ePtr->ce[i][j] = ePtr->elementProps.damping * meBuildingBlock[i][j];
				}
			}
		}
	}
}

