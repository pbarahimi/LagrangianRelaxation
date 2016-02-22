import java.util.ArrayList;
import java.util.List;

import gurobi.*;
import Jama.Matrix;


public class BBNode implements Comparable<BBNode> {
	public int ind;
	public int noOfFixedVars;
	public Matrix Us;
	public GRBVar varToFix = null;
	public List<GRBVar> varsFixedTo1 = new ArrayList<GRBVar>();
	public List<GRBVar> varsFixedTo0 = new ArrayList<GRBVar>();
	public boolean fathom = false;
	public LB lb = new LB();
	public double UB = -1 * GRB.INFINITY;
	public double bestLB;
	public List<String> selectedLocations = new ArrayList<String>();
	public double gap = GRB.INFINITY;
//	public double optGap = 0.01;
	
	private Matrix nodesBestUs;
	private double nodeBestUB;
	private LB nodesBestLB = new LB(-1*GRB.INFINITY);
	public List<String> nodesBestSelectedLocations = new ArrayList<String>();

	
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
	public BBNode(int ind, int N, Matrix Us, LB LB, double UB, List<String> selectedLocations) throws GRBException{
		this.ind = ind;
		this.Us = Us;
		this.lb = LB;
		this.UB = UB; 
		this.selectedLocations = selectedLocations;
		updateGap();
	}
	
	/**
	 * Constructor
	 * @param ind
	 * @param N
	 * @param Us
	 * @param Ds
	 * @param ds
	 * @param OP
	 * @param SP
	 * @param bestLB
	 * @param fixTo1
	 * @param fixTo0
	 * @param varToFix
	 * @param valueToFix
	 * @throws GRBException
	 */
	public BBNode(int ind, int N, Matrix Us, Matrix Ds, Matrix ds, GRBModel OP, GRBModel SP, double bestLB,
			List<GRBVar> fixTo1, List<GRBVar> fixTo0, GRBVar varToFix, boolean valueToFix) throws GRBException{
		this.ind = ind;
		this.Us = Us;
		this.bestLB = bestLB;
		this.varsFixedTo0.addAll(fixTo0);
		this.varsFixedTo1.addAll(fixTo1);
		
		if(valueToFix)
			this.varsFixedTo1.add(varToFix);
		else
			this.varsFixedTo0.add(varToFix);
		
		LR_Main.fixVar(varsFixedTo0, false, SP);
		LR_Main.fixVar(varsFixedTo1, true, SP);
		double miu;
		int cntr = 0;
		SP.optimize();
		if (SP.get(GRB.IntAttr.Status) != 3){  // Model feasibility check
			UB = SP.get(GRB.DoubleAttr.ObjVal);
//			System.out.println("miu: " + miu + " - Itr" + k + ": LB = " + lb.value + " - UB = " + UB);

			while(/*this.gap > optGap && */cntr<5){
//				miu = LR_Main.updateMiu(miu, k, N);
				miu = LR_Main.updateMiuC(Ds, ds);
				this.Us = LR_Main.updateU(SP, miu, this.Us, ds, Ds);		// Update Lagrangian multipliers
				LR_Main.updateObjCoeffs(SP, this.Us, ds, Ds);	// Update obj fun coefficients
				SP.optimize();
				UB = SP.get(GRB.DoubleAttr.ObjVal);
				this.selectedLocations = new ArrayList<String>();
				for (int i = 1 ; i <= N ; i++){
					GRBVar var = SP.getVar(SP.get(GRB.IntAttr.NumVars)-i);
					if (var.get(GRB.DoubleAttr.X) > 0)
						this.selectedLocations.add(var.get(GRB.StringAttr.VarName));					
				}
				this.lb = LR_Main.obtainLB(OP, SP);
				LR_Main.updateLB(lb.value);
				cntr++;							// Counter update
				updateGap();	
				updateBestResult(this.lb, this.UB, this.Us, this.selectedLocations);
				System.out.println("miu: " + miu + " - Itr" + cntr + ": LB=" + lb.value + " - UB= " + UB + " - gap= " + this.gap);
			}
			updateNode(nodeBestUB, nodesBestLB, nodesBestUs, nodesBestSelectedLocations);
		}else{
			this.fathom = true;   // If the node is infeasible, it's fathomed.
		}			
		noOfFixedVars = this.varsFixedTo0.size() + this.varsFixedTo1.size();		
	}
	
	/**
	 * Updates the node attributes to the best found LB, UB, Lagrangian multipliers and Locations.
	 * 
	 * @param bestBB_UB
	 * @param bestBB_LB
	 * @param bestBB_Us
	 * @param bestBB_selectedLocations
	 */
	private void updateNode(double bestBB_UB, LB bestBB_LB, Matrix bestBB_Us, List<String> bestBB_selectedLocations){
		if (this.lb.value <= bestBB_LB.value){
			this.lb = bestBB_LB;
			this.UB = bestBB_UB;
			this.Us = bestBB_Us;
			this.selectedLocations = bestBB_selectedLocations;
		}
	}
	
	/**
	 * Checks the current iteration and updates the best lower bound, upper bound, and locations found so far, if the 
	 * LB of the current iteration is better than the best LB found so far.
	 * 
	 * @param lb
	 * @param ub
	 * @param Us
	 * @param selectedLocations
	 */
	private void updateBestResult(LB lb, double ub, Matrix Us, List<String> selectedLocations){
		if (lb.value >= this.nodesBestLB.value) {
			this.nodesBestLB = lb;
			this.nodeBestUB = ub;
			this.nodesBestUs = Us;
			this.nodesBestSelectedLocations = new ArrayList<String>();
			this.nodesBestSelectedLocations.addAll(selectedLocations);
		}
	}
	
	/**
	 * updates the gap between the UB and LB found in the node.
	 */
	public void updateGap(){
		this.gap = Math.abs((UB - lb.value)/lb.value);
	}
	
	/**
	 * Checks the upper bound of the node with the given LB and fathoms the node
	 * if the UB is less than or equal to LB
	 * @param LB - Lower Bound
	 */
	public void updateFathom(int P, int N){
		if (this.lb.value<bestLB 
				|| this.varsFixedTo1.size() >= P
				|| (N-this.varsFixedTo0.size()) < P) 
			this.fathom = true;
	}
	
	@Override
	public int compareTo(BBNode other){
		double temp = this.lb.value - other.lb.value;
		if (temp>0) return 1;
		else if (temp<0) return -1;
		else return 0;
	}
}
