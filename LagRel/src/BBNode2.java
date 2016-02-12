import java.util.ArrayList;
import java.util.List;

import gurobi.*;
import Jama.Matrix;


public class BBNode2 {
	public int ind;
	public int noOfFixedVars;
	public GRBVar varToFix;
	public List<GRBVar> varsFixedTo1 = new ArrayList<GRBVar>();
	public List<GRBVar> varsFixedTo0 = new ArrayList<GRBVar>();
	public boolean fathom = false;
	public double objVal;
	public List<String> selectedLocations = new ArrayList<String>();
	
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
	public BBNode2(int ind, double bestObj) throws GRBException{
		this.ind = ind;
		this.objVal = bestObj;
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
	public BBNode2(int ind, int N, int nVar, GRBModel OP, double bestObj,
			List<GRBVar> fixTo1, List<GRBVar> fixTo0, GRBVar varToFix, boolean valueToFix) throws GRBException{
		this.ind = ind;
		this.varsFixedTo0.addAll(fixTo0);
		this.varsFixedTo1.addAll(fixTo1);
		
		if(valueToFix)
			this.varsFixedTo1.add(varToFix);
		else
			this.varsFixedTo0.add(varToFix);
		
		LR_Main.fixVar(varsFixedTo0, false, OP);
		LR_Main.fixVar(varsFixedTo1, true, OP);
		OP.optimize();
		if (OP.get(GRB.IntAttr.Status) != 3){  // Model feasibility check
			this.objVal = OP.get(GRB.DoubleAttr.ObjVal);
			for (int i = 1 ; i <= N ; i++){
				GRBVar var = OP.getVar(nVar-i);
				if (var.get(GRB.DoubleAttr.X) > 0)
					this.selectedLocations.add(var.get(GRB.StringAttr.VarName));					
			}
		}else{
			this.fathom = true;
		}			
//		System.out.println("# of constraints: " + SP.get(GRB.IntAttr.NumConstrs));
		noOfFixedVars = this.varsFixedTo0.size() + this.varsFixedTo1.size();		
	}
	
	/**
	 * Checks the upper bound of the node with the given LB and fathoms the node
	 * if the UB is less than or equal to LB
	 * @param LB - Lower Bound
	 */
	public void updateFathom(double bestObj, int P, int N){
		if (this.objVal < bestObj ||this.varsFixedTo1.size() >= P|| (N-this.varsFixedTo0.size()) < P) 
			this.fathom = true;
	}
}
