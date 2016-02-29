import java.util.ArrayList;
import java.util.List;

import gurobi.*;
import Jama.Matrix;


public class BBNode2 implements Comparable<BBNode2> {
	public int ind;
	public int noOfFixedVars;
	public List<Double> Us;
	public List<GRBVar> slacks;
	public GRBVar varToFix;
	public List<GRBVar> varsFixedTo1 = new ArrayList<GRBVar>();
	public List<GRBVar> varsFixedTo0 = new ArrayList<GRBVar>();
	public boolean fathom = false;
	public LB lb = new LB(-1*GRB.INFINITY);
	public double UB = GRB.INFINITY;
	public double bestLB;
	public double bestUB = GRB.INFINITY;
	public List<String> selectedLocations = new ArrayList<String>();
	public double gap = GRB.INFINITY;
	public double optGap = 0.01;
	public double epsilon;
	
	private List<Double> nodesBestUs;
	private double nodeBestUB;
	private LB nodesBestLB = new LB(-1*GRB.INFINITY);
	public List<String> nodesBestSelectedLocations = new ArrayList<String>();

	
	/**
	 * Constructor
	 * @param ind
	 * @param slacks
	 * @param Us
	 * @param LB
	 * @param UB
	 * @param selectedLocations
	 * @throws GRBException
	 */
	public BBNode2(int ind, List<GRBVar> slacks, List<Double> Us, LB LB, double UB, List<String> selectedLocations, double epsilon) throws GRBException{
		this.ind = ind;
		this.slacks = slacks;
		this.Us = Us;
		this.nodesBestUs = Us;
		this.bestLB = LB.value;
		this.UB = UB; 
		this.selectedLocations = selectedLocations;
		this.epsilon = epsilon;
		updateGap();
	}
	
	/**
	 * Constructor
	 * @param ind
	 * @param N
	 * @param slacks
	 * @param Us
	 * @param OP
	 * @param SP
	 * @param bestLB
	 * @param fixTo1
	 * @param fixTo0
	 * @param varToFix
	 * @param valueToFix
	 * @throws GRBException
	 */
	public BBNode2(int ind, int N, List<GRBVar> slacks, List<Double> Us, GRBModel OP, GRBModel SP, double bestLB,
			List<GRBVar> fixTo1, List<GRBVar> fixTo0, GRBVar varToFix, boolean valueToFix, double epsilon, int itrNum) throws GRBException{
		this.ind = ind;
		this.slacks = slacks;
		this.Us = Us;
		this.nodesBestUs = this.Us;
		this.bestLB = bestLB;
		this.epsilon = epsilon;
		this.varsFixedTo0.addAll(fixTo0);
		this.varsFixedTo1.addAll(fixTo1);
		if(valueToFix)
			this.varsFixedTo1.add(varToFix);
		else
			this.varsFixedTo0.add(varToFix);
		
		LR_Main2.fixVar(varsFixedTo0, false, SP);
		LR_Main2.fixVar(varsFixedTo1, true, SP);
		double miu;
		int cntr = 0;
		int k = 0;
		SP.optimize();
		if (SP.get(GRB.IntAttr.Status) != 3){  // Model feasibility check
			UB = SP.get(GRB.DoubleAttr.ObjVal);
			while(/*this.gap > optGap &&*/ k<itrNum){
				if (LR_Main2.dissectEpsilon(cntr, 3)){
					cntr=0;
					this.epsilon = this.epsilon/2;
				}
				miu = LR_Main2.updateMiuC(this.slacks, this.epsilon);
				LR_Main2.updateU(this.Us, this.slacks, miu);		// Update Lagrangian multipliers
				LR_Main2.updateObjCoeffs(this.slacks, this.Us);	// Update obj fun coefficients
				SP.optimize();
				UB = SP.get(GRB.DoubleAttr.ObjVal);
				if (UB<this.bestUB){
					this.bestUB = UB;
					cntr = 0;
				}else{
					cntr++;
				}
				this.selectedLocations = new ArrayList<String>();
				for (int i = 1 ; i <= N ; i++){
					GRBVar var = SP.getVar(SP.get(GRB.IntAttr.NumVars)-i);
					if (var.get(GRB.DoubleAttr.X) > 0)
						this.selectedLocations.add(var.get(GRB.StringAttr.VarName));					
				}
				lb = LR_Main2.obtainLB(OP, SP);
				LR_Main2.updateLB(lb.value);
				k++;							// Counter update
				updateGap();	
				updateBestResult(this.lb, this.UB, this.Us, this.selectedLocations);
				System.out.println("miu: " + miu + " - Itr" + k + ": LB=" + lb.value + " - UB= " + UB + " - gap= " + this.gap +  " - eps= " + LR_Main.epsilon);
			}
			updateNode(nodeBestUB, nodesBestLB, nodesBestUs, nodesBestSelectedLocations);
		}else{
			this.fathom = true;
		}			

		noOfFixedVars = this.varsFixedTo0.size() + this.varsFixedTo1.size();		
	}
	
	private void updateNode(double bestBB_UB, LB bestBB_LB, List<Double> bestBB_Us, List<String> bestBB_selectedLocations){
		if (this.lb.value <= bestBB_LB.value){
			this.lb = bestBB_LB;
			this.UB = bestBB_UB;
			this.Us = bestBB_Us;
			this.selectedLocations = bestBB_selectedLocations;
		}
	}
	
	private void updateBestResult(LB lb, double ub, List<Double> Us, List<String> selectedLocations){
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
		this.gap = Math.abs((UB - lb.value)/bestLB);
	}
	
	/**
	 * Checks the upper bound of the node with the given LB and fathoms the node
	 * if the UB is less than or equal to LB
	 * @param LB - Lower Bound
	 */
	public void updateFathom(int P, int N){
		if (this.lb.value < bestLB 
				|| this.varsFixedTo1.size() >= P
				|| (N-this.varsFixedTo0.size()) < P) this.fathom = true;
	}
	
	@Override
	public int compareTo(BBNode2 other){
		double temp = this.lb.value - other.lb.value;
		if (temp>0) return -1;
		else if (temp<0) return 1;
		else return 0;
	}
}