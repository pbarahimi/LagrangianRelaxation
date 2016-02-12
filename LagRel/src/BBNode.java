import java.util.ArrayList;
import java.util.List;

import gurobi.*;
import Jama.Matrix;


public class BBNode {
	public int ind;
	public int noOfFixedVars;
	public Matrix Us;
	public GRBVar varToFix;
	public List<GRBVar> varsFixedTo1 = new ArrayList<GRBVar>();
	public List<GRBVar> varsFixedTo0 = new ArrayList<GRBVar>();
	public boolean fathom = false;
	public LB lb = new LB();
	public double UB = -1 * GRB.INFINITY;
	public List<String> selectedLocations = new ArrayList<String>();
	public double gap = 0;
	public double globalGap = 0;
	
	/**
	 * Constructor 
	 * @param ind - node index
	 * @param N - number of nodes
	 * @param Us - Lagrangian Multipliers
	 * @param LB - Lower Bound
	 * @param UB - Upper Bound
	 * @param OP - 
	 * @throws GRBException 
	 */
	public BBNode(int ind, int N, Matrix Us, LB LB, double UB) throws GRBException{
		this.ind = ind;
		this.Us = Us;
		this.lb = LB;
		this.UB = UB; 
	}
	
	/**
	 * Constructor
	 * @param ind - node index
	 * @param N - number of nodes
	 * @param Us - matrix of Lagrangian multipliers
	 * @param Ds - matrix of coefficients
	 * @param ds - right hand side matrix
	 * @param OP - original problem
	 * @param SP - sub-problem
	 * @param fixTo1 - list of variables to be set to 1
	 * @param fixTo0 - list of variables to be set to 0
	 * @throws GRBException
	 */
	public BBNode(int ind, int N, Matrix Us, Matrix Ds, Matrix ds, GRBModel OP, GRBModel SP, double lowerBound,
			List<GRBVar> fixTo1, List<GRBVar> fixTo0, GRBVar varToFix, boolean valueToFix) throws GRBException{
		this.ind = ind;
		this.Us = Us;
		this.lb.value = lowerBound;
		this.varsFixedTo0.addAll(fixTo0);
		this.varsFixedTo1.addAll(fixTo1);
		
		if(valueToFix)
			this.varsFixedTo1.add(varToFix);
		else
			this.varsFixedTo0.add(varToFix);
		
		LR_Main.fixVar(varsFixedTo0, false, SP);
		LR_Main.fixVar(varsFixedTo1, true, SP);
		double miu = 1;
		int k = 0;
		SP.optimize();
		if (SP.get(GRB.IntAttr.Status) != 3){  // Model feasibility check
			UB = SP.get(GRB.DoubleAttr.ObjVal);
//			System.out.println("miu: " + miu + " - Itr" + k + ": LB = " + lb.value + " - UB = " + UB);

			double terCri=1; 	// Termination criteria
			while(Math.abs(UB-lb.value)>lb.value*0.1 && terCri!=0 && miu>0.01){
				miu = LR_Main.updateMiu(miu, k, N);
				this.Us = LR_Main.updateU(SP, miu, this.Us, ds, Ds);		// Update Lagrangian multipliers
				LR_Main.updateObjCoeffs(SP, this.Us, ds, Ds);	// Update obj fun coefficients
				SP.optimize();
//				terCri = UB - SP.get(GRB.DoubleAttr.ObjVal);
				UB = SP.get(GRB.DoubleAttr.ObjVal);
//				LR_Main.printSol(SP);
				this.selectedLocations.clear();
				for (int i = 1 ; i <= N ; i++){
					GRBVar var = SP.getVar(SP.get(GRB.IntAttr.NumVars)-i);
					if (var.get(GRB.DoubleAttr.X) > 0)
						this.selectedLocations.add(var.get(GRB.StringAttr.VarName));					
				}
				lb = LR_Main.obtainLB(OP, SP);
				k++;							// Counter update
//				System.out.println("miu: " + miu + " - Itr" + k + ": LB=" + lb.value + " - UB= " + UB);
			}
			updateGap();
		}else{
			this.fathom = true;
		}			

//		System.out.println("# of constraints: " + SP.get(GRB.IntAttr.NumConstrs));
		noOfFixedVars = this.varsFixedTo0.size() + this.varsFixedTo1.size();		
	}
	
	/**
	 * updates the gap between the UB and LB found in the node.
	 */
	public void updateGap(){
		this.gap = UB - lb.value; 
	}
	
	/**
	 * Checks the upper bound of the node with the given LB and fathoms the node
	 * if the UB is less than or equal to LB
	 * @param LB - Lower Bound
	 */
	public void updateFathom(double bestGlobalGap, double LB, int P, int N){
		if (this.UB<=LB || this.gap > bestGlobalGap) this.fathom = true;
		else if (this.varsFixedTo1.size() >= P
				|| (N-this.varsFixedTo0.size()) < P) this.fathom = true;
	}
}
