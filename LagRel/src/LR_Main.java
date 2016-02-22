import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import Jama.Matrix;
import gurobi.*;

public class LR_Main {
	// formulation attributes
	private static double alpha = 0.2;
	private static double[][] tmpFlows = MyArray.read("w.txt");
	private static double[][] coordinates = MyArray.read("coordinates.txt");
	private static double[][] distances = Distance.get(coordinates);
	private static int N = tmpFlows.length;
	private static double[][] flows = new double[N][N];
	private static int D = 2; // Maximum number of simultaneous disruptions
	private static int R = (int) Math.pow(2, D + 1) - 2; // Largest index in the
															// full binary tree
	private static double q = 0.3;
	private static int P = 3; // number of hubs to be located
	private static int M = N * R; // the big M
	private static int nVar = (int) (N + (Math.pow(N, 3)*(N-1)*(R+1)/2)); // Total number of variables in the model
	
	// Lagrangian Relaxation algorithm attributes
	private static GRBModel SP;	// Sub problem
	private static GRBModel OP;	// Original problem
	private static Matrix Us;
	private static Matrix Ds;
	private static Matrix ds;
	private static Matrix C = new Matrix(1,(int) ((Math.pow(N, 3)*(N-1)*(R+1)/2)+N));
	private static int k = 0;
	private static double miu = 6;
	private static double epsilon = 2;
	private static double UB;
	private static LB LB;
	private static double bestUB = GRB.INFINITY;
	private static double bestLB = -1 * GRB.INFINITY;
	private final static double gap = 0.05; // Termination criteria: the difference between the 
											// last obtained objective function value and the second to the last
	
	// Branch and Bound attributes
//	private static List<BBNode> exploredNodes = new ArrayList<BBNode>();
	private static List<BBNode> unexploredNodes = new ArrayList<BBNode>();
	private static BBNode bestNode;
	private static PriorityQueue<BBNode> BBNodeList = new PriorityQueue<BBNode>();
	
	// Variable fixing attributes
	private static double z; // incumbent value for variable fixing
	private static List<GRBVar> N0;
	private static List<GRBVar> N1;
	private static List<GRBVar> fixTo1;
	private static List<GRBVar> fixTo0;
	
	
	/*
	 * Formulation methods
	 */
	
	/**
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @param m
	 * @return operating probability of a route
	 */
	public static double Q(int i, int k, int m, int j) {
		double result = q;
		if (k!=i && j!=m)
			result = q+q(k,m);
		else if (m!=j)
			result = q(i,m);
		else if (k!=i)
			result = q(j,k);
		else if (i==k && j==m)
			result = 0;
		else 
			System.out.println("Not include in the Q(i,k,m,j)!");
		return result;
	}

	/**
	 * Cikmj
	 * 
	 */
	private static double Cikmj(int i, int k, int m, int j) {
		double cost = distances[i][k] + (1 - alpha) * distances[k][m]
				+ distances[m][j];
		/*
		 * double cost = collCost * distances[i][k] + transCost *
		 * distances[k][m] + distCost * distances[m][j];
		 */
		return cost;
	}

	/**
	 * q(k,m)
	 */
	private static double q(int k, int m) {
		if (k == m)
			return 0;
		else
			return q;
	}
	
	/*
	 * Lagrangian relaxation methods
	 *
	 */
	
	/**
	 * Updates the bestLB found so far.
	 * @param x
	 */
	protected static void updateLB(double x){
		if (x > bestLB) bestLB = x;
	}
	
	/**
	 * Compares the new UB found with the UB found so far and replaces it 
	 * if its better.
	 * @param UB
	 */
	protected static int updateUB(int counter, double UB){
		if (UB<bestUB){
			bestUB = UB;
			counter = 0;
		}else{
			counter++;
		}
		return counter;
	}
	
	/**
	 * Halves the epsilon and resets the counter if the threshold is met
	 */
	protected static int dissectEpsilon(int counter, int threshold){
		if (counter >= threshold){
			epsilon = epsilon/2;
			counter = 0;
		}
		return counter;		 
	}
	
	/**
	 * prints the solution
	 * @param model
	 * @throws GRBException
	 */
	protected static void printSol(GRBModel model) throws GRBException{
		for (GRBVar var : model.getVars()) 
			if (var.get(GRB.DoubleAttr.X)>0 /* && var.get(GRB.StringAttr.VarName).contains("y")*/) 
				System.out.println(var.get(GRB.StringAttr.VarName) + " : " + var.get(GRB.DoubleAttr.X));
		System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));
	}
	
	/**
	 * prints the solution
	 * @param model
	 * @throws GRBException
	 */
	protected static String printSol2(GRBModel model) throws GRBException{
		String output = "values: ";
		for (GRBVar var : model.getVars()) 
			if (var.get(GRB.DoubleAttr.X)>0  && var.get(GRB.StringAttr.VarName).contains("y")){ 
				output = output.concat(var.get(GRB.StringAttr.VarName)+ ",");
			}
		return output;
	}
	
	/**
	 * Updates the objective function coefficients based on the Lagrangian Multipliers put in.
	 * @param U
	 * @param d
	 * @param D
	 * @throws GRBException
	 */
	protected static void updateObjCoeffs(GRBModel SP, Matrix U, Matrix d, Matrix D) throws GRBException{
		double[][] temp =  C.getArrayCopy();
		Matrix updatedC = new Matrix(temp);
		double constant = 0;
//		for (int i = 0 ; i<U.size() ; i++){
			updatedC = updatedC.minus(U.transpose().times(D));			
			constant += (U.transpose().times(d)).get(0, 0);
//		}
		GRBVar var;
		for (int i = 0 ; i < nVar - N ; i++){
			var = SP.getVar(i);
			var.set(GRB.DoubleAttr.Obj, updatedC.getArray()[0][getCol(var)]);
		}
		for (int i = 0 ; i < N ; i++){
			var = SP.getVar(nVar-i-1);
			var.set(GRB.DoubleAttr.Obj, updatedC.getArray()[0][getLocCol(var)]);
		}
		SP.set(GRB.DoubleAttr.ObjCon, constant);		
	}

	/**
	 * updates the step-size based on the iteration number
	 * @param miu
	 * @param k
	 * @return
	 */
	protected static double updateMiu(double miu, int k, int N){
		int itr = (int) Math.floor(k);
		miu = miu* Math.pow(0.99,itr);
		return miu;
	}
	
	/**
	 * Concatenates a list of matrices horizontally
	 * @param x
	 * @return
	 */
	private static Matrix concatH(List<Matrix> x){
		int colDim = x.get(0).getColumnDimension(); // column dimension of the output matrix
		int rowDim = 0;	// row dimension of the output matrix
		for (Matrix i:x){
			rowDim += i.getRowDimension();
		}
		Matrix output = new Matrix(rowDim, colDim);
		int i0;
		int i1 = -1;
		for (Matrix i:x){
			i0 = i1+1;
			i1 = i1 + i.getRowDimension();
			output.setMatrix(i0, i1, 0, colDim-1, i);			
		}
		return output;
	}
	
	/**
	 * Gets the column number of routing variables.
	 * @param var
	 * @return int
	 * @throws GRBException
	 */
	private static int getCol(GRBVar var) throws GRBException{
		VarIndices varIndices = new VarIndices(var);
		int _R = R+1;
		int lambda = N * _R * (N-1) + _R*(N-1) + _R;    // Index j increments every lambda iterations.
		int output = 0;
		for (int i = 0 ; i < varIndices.i; i++){
			output += (N-i-1) * lambda;
		}		
		output += (varIndices.j - varIndices.i -1) * lambda + (N*_R*varIndices.k) + (_R*varIndices.m) + varIndices.r;		
		return output;
	}
	
	/**
	 * Get the column number of location variables.
	 * @param var
	 * @return int
	 * @throws NumberFormatException
	 * @throws GRBException
	 */
	private static int getLocCol(GRBVar var) throws NumberFormatException, GRBException{
		int output = (N*(N-1)/2*N*N*(R+1));
		return output + Integer.parseInt(var.get(GRB.StringAttr.VarName));
	}
	
	/**
	 * Updates the miu based on the  rule (c)
	 * @param model
	 * @param varHash
	 * @param UB
	 * @param currentU
	 * @param epsilon
	 * @return
	 * @throws GRBException 
	 */
	protected static double updateMiuC(Matrix D, Matrix d) throws GRBException{
//		double gap1 = delta2 - delta1;
		Matrix X = getVarsMatrix(SP);
		Matrix denuminatorMat = D.times(X.transpose());		
		denuminatorMat = d.minus(denuminatorMat);
		double denuminator = Math.pow(denuminatorMat.normF(),2);	 
		double output = epsilon*(SP.get(GRB.DoubleAttr.ObjVal) - LB.value ) / denuminator;
//		delta1 = delta2;
//		delta2 = SP.get(GRB.DoubleAttr.ObjVal) - UB;
//		double gap2 = delta2 - delta1;
//		if (gap2>gap1) UB += Math.abs(0.05 * UB);
		return output;
	}
	
	/**
	 * Updates the Lagrangian Multipliers according to the miu.
	 * @param U
	 * @param d
	 * @param D
	 * @param miu
	 * @param model
	 * @param varHash
	 * @return
	 * @throws GRBException
	 */
	protected static Matrix updateU(GRBModel SP, double miu, Matrix U, Matrix d, Matrix D) throws GRBException{
		Matrix X = getVarsMatrix(SP);
		Matrix output = D.times(X.transpose());		
		output = d.minus(output);
		output = U.minus(output.times(miu));
		double[][] u_k = output.getArray();
		for (int i=0;i<u_k.length;i++){
			for (int j=0;j<u_k[0].length;j++){
				u_k[i][j] = (u_k[i][j] > 0) ? u_k[i][j]:0;
			}
		}
		Matrix result = new Matrix(u_k);
		return result;
	}
	
	/**
	 * converts the array of Gurobi variables into a 1xn matrix.
	 * @param vars
	 * @return
	 * @throws GRBException
	 */
	private static Matrix getVarsMatrix(GRBModel model) throws GRBException{
		Matrix output = new Matrix(1, model.get(GRB.IntAttr.NumVars));
		for (int i = 0 ; i < nVar-N ; i++){
			GRBVar var = model.getVar(i);
			output.set(0, getCol(var), var.get(GRB.DoubleAttr.X));			
		}
		for (int i = 0 ; i < N ; i++){
			GRBVar var = model.getVar(nVar - i - 1);
			output.set(0, getLocCol(var), var.get(GRB.DoubleAttr.X));
		}
		return output;
	}
	
	/**
	 * Checks for the termination criteria. If the difference between the last objValue
	 * and the second last is less than the "gap" the algorithm stops.
	 * @param newObj
	 * @param oldObj
	 * @return true if the termination criteria is not met.
	 */
	private static boolean run(double newObj, double oldObj){
		if (Math.abs(newObj-oldObj) > gap) return true;
		else return false;		
	}
	
	/*
	 * Branch and Bound methods
	 */
	
	/**
	 * Calculates the flow that goes through each hub and
	 * select the one that is selected by the OP and has 
	 * the most flow going through.
	 * 
	 * @param x - route variables
	 * @param y - location variables
	 * @return GRBVar
	 * @throws GRBException
	 */
	public static GRBVar getBranchVar(GRBVar[][][][][] x, List<GRBVar> y) throws GRBException{
		/*
		 * Calculating the flow that goes through each hub and
		 * selecting the one that is selected by the OP and has 
		 * the most flow going through it.
		 */
		double[] flowSumOfHubs = new double[N];
		
		for (int r = 0; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int k = 0; k < N; k++) {
						for (int m = 0 ; m < N ; m++){
							flowSumOfHubs[k] += x[i][k][m][j][r].get(GRB.DoubleAttr.Obj);
							flowSumOfHubs[m] += x[i][k][m][j][r].get(GRB.DoubleAttr.Obj);
						}
					}
				}
			}
		}
		
		double largestFlow = -1 * GRB.INFINITY;
		GRBVar branchVar = y.get(0);
		for ( GRBVar var : y ){		
			if (var.get(GRB.DoubleAttr.X)==1 && flowSumOfHubs[Integer.parseInt(var.get(GRB.StringAttr.VarName))]>largestFlow){
				branchVar = var;
			}
		}
		return branchVar;
	}
	
	/**
	 * Updates the LB by comparing the LB obtained at node
	 * @param node
	 */
	public static void updateLB(BBNode node){
		if (node.lb.value>LB.value)
			LB.value = node.lb.value;
	}
	
	/**
	 * Takes the location variables and a node, removes the fixed variables
	 * from the list and returns the unfixed location variables
	 * @param y -  array of location variables
	 * @param node - B&B node
	 * @return
	 */
	public static List<GRBVar> getUnfixedVars(GRBVar[] y, BBNode node){
		List<GRBVar> output = new ArrayList<GRBVar>();
		for (int i=0; i<y.length ; i++)
			output.add(y[i]);
		
		output.removeAll(node.varsFixedTo0);
		output.removeAll(node.varsFixedTo1);
		return output;		
	}
	
	/**
	 * Gets a node and replaces it as the best node if
	 * its upper bound is greater than the current best node upper bound.
	 * @param node
	 */
	public static void updateBestNode(BBNode node){
		if (node.UB > bestNode.UB)
			bestNode = node;
	}
	
	/**
	 * Prints B&B node attributes
	 * @param node
	 * @throws GRBException
	 */
	public static void printNode(BBNode node) throws GRBException{
		System.out.println("fathomed: " + node.fathom);
		System.out.println("node UB: " + node.UB);
		System.out.println("node LB: " + node.lb.value);
		System.out.println("Gap: " + node.gap);
		if (node.varToFix != null) System.out.println("var to be fixed: " + node.varToFix.get(GRB.StringAttr.VarName));
		System.out.println("Locations selected: " + node.selectedLocations);
		System.out.println("-----------------------------");
	}
	
	/*
	 * Variable fixing methods
	 */
	
	/**
	 * Calculates the benefits of locating a hub in any of the nodes and divides the location variables
	 * into selected and unselected location and put them in min and max heap respectively.
	 *  
	 * @param model
	 * @param x - routing variables
	 * @param y - location variables
	 * @return
	 * @throws GRBException
	 */
	public static BenefitList computeBenefits(GRBModel model, GRBVar[][][][][] x, GRBVar[] y) throws GRBException{
		BenefitList output = new BenefitList();
		double benefit;
		
		// computing benefits
		double temp;
		double value;
		for (int k = 0 ; k < N ; k++){
			benefit = 0;
			
			for (int i = 0 ; i < N ; i++){
				for (int j = i+1 ; j < N ; j++){
					
					value = -1*GRB.INFINITY;
					for (int m = 0 ; m < N ; m++){
						for (int r = 0 ; r < R ; r++){
							temp = Math.max(x[i][k][m][j][r].get(GRB.DoubleAttr.Obj), x[i][m][k][j][r].get(GRB.DoubleAttr.Obj));
							value = (value > temp) ? value : temp;							
						}
					}
					benefit += Math.max(0,value);					
				}
			}
			output.add(y[k], k, benefit);
		}
		
		return output;
	}
	
	public static void fixVar2(BenefitList list, GRBVar[] y, double UB, double LB){
		for (BenefitList.Benefit location : list.selectedLocs)
			if (UB - location.benefit + list.unselectedLocs.peek().benefit < LB) fixTo1.add(y[location.index]); 
		
		for (BenefitList.Benefit location : list.unselectedLocs)
			if (UB + location.benefit - list.selectedLocs.peek().benefit < LB) fixTo0.add(y[location.index]);		
	}
	
	/**
	 * Adds a constraint to the model that fixes the variable and assigns a name to the constraint 
	 * @param var
	 * @param i
	 * @param s
	 * @throws GRBException
	 */
	public static void fixVar(GRBVar var, double i, String s) throws GRBException{
		GRBLinExpr expr = new GRBLinExpr();
		expr.addTerm(1, var);
		OP.addConstr(expr, GRB.EQUAL, i, s);
	}
	
	/**
	 * Returns the sum of elements in a matrix.
	 * @param m
	 * @return
	 */
	private static double sumElements(Matrix m){
		double[][] arr = m.getArray();
		double output = 0;
		for (int i = 0 ; i < arr.length ; i++){
			for (int j = 0 ; j < arr[0].length ; j++){
				output += arr[i][j];
			}
		}
		return output;
	}
	
	/**
	 * Takes a set of GRBVar and assigns them to N1 and N0.
	 * @param vars
	 * @throws GRBException
	 */
	private static void setN1N0(GRBVar[] vars) throws GRBException{
		N0 = new ArrayList<GRBVar>();
		N1 = new ArrayList<GRBVar>();
		int col;
		double criteria;
		for (GRBVar var:vars){
			col = getCol(var); // column index of the variable
			criteria = C.get(0, col)-sumElements(Us.arrayTimes(Ds.getMatrix(0, Ds.getRowDimension()-1, col, col)));
			System.out.println(criteria);
			if (criteria > 0)
				N1.add(var);
			else if (criteria < 0)
				N0.add(var);
		}
	}
	
	/**
	 * Defines variables to be fixed
	 * @param N1
	 * @param N0
	 * @param z The incumbent value
	 * @throws GRBException 
	 */
	private static void varsTobeFixed(List<GRBVar> N1, List<GRBVar> N0, double z) throws GRBException{
		fixTo0 = new ArrayList<GRBVar>();
		fixTo1 = new ArrayList<GRBVar>();
		int rowSize = Us.getRowDimension()-1; // Number of rows in Us, ds and Ds.
		double criteria;
		double statement1 = sumElements(Us);  // The first statement in the variable fixing formula
		double statement2 = 0;  // The second statement in the variable fixing formula
		for (GRBVar var : N0){
			int col = getCol(var);
			statement2 += var.get(GRB.DoubleAttr.Obj) /*C.get(0,col)*/ - sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));			
		}
		
		// Defining variables to be fixed to 0
		for (GRBVar var : N1){			
			int col = getCol(var);
			criteria = statement1 + statement2 + var.get(GRB.DoubleAttr.Obj) /*C.get(0, col)*/ - sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));
			System.out.println("cri: "+ criteria);
			if (criteria >= z) fixTo0.add(var);
			}
		// Defining variables to be fixed to 1
		double statement2_2;
		for (GRBVar var : N0){			
			int col = getCol(var);
			statement2_2 = statement2 - C.get(0,col) + sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));
			criteria = statement1 + statement2_2;	
			System.out.println("cri: "+ criteria);
			if (criteria >= z) fixTo1.add(var);
		}		
	}
	
	/**
	 * Relaxes the binary variables in the model.
	 * @param model
	 * @throws GRBException
	 */
	private static void relaxModel(GRBModel model) throws GRBException{
		for (GRBVar var : model.getVars())
			var.set(GRB.CharAttr.VType, 'C');		
	}
	
	/**
	 * Sets all the variables in the model to binaries
	 * @param model
	 * @throws GRBException
	 */
	private static void setBinaries(GRBModel model) throws GRBException{
		for (GRBVar var : model.getVars())
			var.set(GRB.CharAttr.VType, 'B');
		model.update();
	}
	
	/**
	 * Obtains lower bound by fixing location variables in the original problem.
	 * @param y
	 * @param y0
	 * @return Lower Bound (double)
	 * @throws GRBException
	 */
	protected static LB obtainLB(GRBModel OP, GRBModel SP) throws GRBException{
		LB lb = new LB();
		int nVar = OP.get(GRB.IntAttr.NumVars);
		for (int i = 1 ; i <= N ; i++)
			fixVar(OP.getVar(nVar-i), SP.getVar(nVar-i).get(GRB.DoubleAttr.X), null);
		OP.optimize();
		lb.value = OP.get(GRB.DoubleAttr.ObjVal);	
		
		// store the name of location variables that are set to 1
		for (int i = 1 ; i<=N ; i++)
			if (OP.getVar(nVar-i).get(GRB.DoubleAttr.X)>0)
				lb.vars.add(OP.getVar(nVar-i).get(GRB.StringAttr.VarName));
		
		// Removing fixing constraints
		relaxVar(OP, N); 
		return lb;
	}
		
	/**
	 * Adds constraints to the model so that the variables are fixed to the "value".
	 * @param vars
	 * @param value
	 * @param model
	 * @throws GRBException
	 */
	protected static void fixVar(List<GRBVar> vars, boolean value, GRBModel model) throws GRBException{
		if (vars.equals(null)) return;
		int rhs = 0;
		if (value) rhs = 1;
		for (GRBVar var : vars){
			GRBLinExpr expr = new GRBLinExpr();
			expr.addTerm(1, var);
			model.addConstr(expr, GRB.EQUAL, rhs, null);
		}		
	}
	
	/**
	 * Relaxes the N last constraints added to the model. The N last constraints are supposed to be
	 * the variable fixing constraints
	 * @param model
	 * @param N
	 * @throws GRBException
	 */
	protected static void relaxVar(GRBModel model, int N) throws GRBException{
		int k;
		for (int i = 0 ; i < N ; i++){
			k = model.get(GRB.IntAttr.NumConstrs) - 1;
			model.remove(model.getConstr(k));
			model.update();
		}	
	}
	
	public static void main(String[] arg) throws GRBException, FileNotFoundException{
		//Filling in the flows matrix asymmetrically 
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
		
		// initializing the model and setting environment attributes
		GRBEnv env = new GRBEnv("RpHND_LR.log");
		SP = new GRBModel(env);
		OP = new GRBModel(env);
		SP.getEnv().set(GRB.IntParam.OutputFlag, 0);
		OP.getEnv().set(GRB.IntParam.OutputFlag, 0);
		
		// initializing variables
		GRBVar[][][][][] x0 = new GRBVar[N][N][N][N][R + 1]; // variables to be added to the original model		
		GRBVar[] y0 = new GRBVar[N]; // variables to be added to the original model	
		GRBVar[][][][][] x = new GRBVar[N][N][N][N][R + 1]; // variables to be added to the Lagrangian dual model (Sub-problem)
		GRBVar[] y = new GRBVar[N];// variables to be added to the Lagrangian dual model (Sub-problem)
		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						for (int r = 0; r <= R; r++){
							String varName = i + "_" + k + "_" + m + "_"
									+ j + "_" + r;
							x[i][k][m][j][r] = SP.addVar(0.0, GRB.INFINITY, 0.0,
									GRB.BINARY, varName);
							x0[i][k][m][j][r] = OP.addVar(0.0, GRB.INFINITY, 0.0,
									GRB.BINARY, varName);
						}
					}
				}
			}
		}		

		for (int i = 0; i < N; i++) {
			y[i] = SP.addVar(0, GRB.INFINITY, 0, GRB.BINARY, ""+i);
			y0[i] = OP.addVar(0, GRB.INFINITY, 0, GRB.BINARY, ""+i);
		}
		
		OP.update();
		SP.update();
		
		// obj function declaration
		GRBLinExpr expr0 = new GRBLinExpr();
		GRBLinExpr expr = new GRBLinExpr();
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						double CoEf = -1 * flows[i][j] * Cikmj(i, k, m, j) * (1 - Q(i,k,m,j));
						expr0.addTerm(CoEf, x0[i][k][m][j][0]);
						expr.addTerm(CoEf, x[i][k][m][j][0]);
						C.set(0, getCol(x[i][k][m][j][0]), CoEf);
					}
				}
			}
		}

		for (int r = 1; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							double CoEf = -1 * flows[i][j] * Cikmj(i, k, m, j) * Math.pow(q, Math.floor(Math.log(r+1)/Math.log(2)));
							expr0.addTerm(CoEf, x0[i][k][m][j][r]);
							expr.addTerm(CoEf, x[i][k][m][j][r]);
							C.set(0, getCol(x[i][k][m][j][r]), CoEf);
						}
					}
				}
			}
		}

		OP.setObjective(expr0, GRB.MAXIMIZE);
		SP.setObjective(expr, GRB.MAXIMIZE);
		

		// Adding constraints
		// Constraint 2
		GRBLinExpr con2_0 = new GRBLinExpr();
		GRBLinExpr con2 = new GRBLinExpr();
		for (int i = 0; i < N; i++) {
			con2_0.addTerm(1, y0[i]);
			con2.addTerm(1, y[i]);
		}
		OP.addConstr(con2_0, GRB.EQUAL, P, "u2");
		SP.addConstr(con2, GRB.EQUAL, P, "u2");

		// Constraint 3
		int cntr = 0;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {

				GRBLinExpr con3_0 = new GRBLinExpr();
				GRBLinExpr con3 = new GRBLinExpr();
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						con3_0.addTerm(1, x0[i][k][m][j][0]);
						con3.addTerm(1, x[i][k][m][j][0]);
					}
				}
				OP.addConstr(con3_0, GRB.EQUAL, 1, "u3_" + i + "_" + j);
				SP.addConstr(con3, GRB.EQUAL, 1, "u3_" + i + "_" + j);				
			}
		}

		// Constraint 4
		cntr = 0;
		for (int r = 0; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int k = 0; k < N; k++) {

						GRBLinExpr con4_0 = new GRBLinExpr();
						GRBLinExpr con4 = new GRBLinExpr();
						for (int m = 0; m < N; m++) {
							con4_0.addTerm(1, x0[i][k][m][j][r]);
							con4.addTerm(1, x[i][k][m][j][r]);
						}
						con4_0.addTerm(-1, y0[k]);
						con4.addTerm(-1, y[k]);
						OP.addConstr(con4_0, GRB.LESS_EQUAL, 0, "u4_" + i + "_" + j + "_" + k + "_" + r);
						SP.addConstr(con4, GRB.LESS_EQUAL, 0, "u4_" + i + "_" + j + "_" + k + "_" + r);
					}
				}
			}
		}

		// Constraint 5
		cntr = 0;
		for (int r = 0; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int m = 0; m < N; m++) {

						GRBLinExpr con5_0 = new GRBLinExpr();
						GRBLinExpr con5 = new GRBLinExpr();
						for (int k = 0; k < N; k++) {
							con5_0.addTerm(1, x0[i][k][m][j][r]);
							con5.addTerm(1, x[i][k][m][j][r]);
						}
						con5_0.addTerm(-1, y0[m]);
						con5.addTerm(-1, y[m]);
						OP.addConstr(con5_0, GRB.LESS_EQUAL, 0, "u5_" + i + "_" + j + "_" + m + "_" + r);
						SP.addConstr(con5, GRB.LESS_EQUAL, 0, "u5_" + i + "_" + j + "_" + m + "_" + r);
					}
				}
			}
		}
		
		
		/*
		 *  Building Lagrangian multipliers and corresponding constraints 
		 */
		// Constraint 6
		Matrix D6 = new Matrix(N*(N-1)*(R+1)/2, nVar);
		Matrix d6 = new Matrix(N*(N-1)*(R+1)/2, 1, M);
		Matrix U6 = new Matrix(N*(N-1)*(R+1)/2, 1, 500);		
		D6 = new Matrix(N*(N-1)*(R+1)/2,nVar);
		cntr = 0;		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= R; r++) {

					GRBLinExpr con6_0 = new GRBLinExpr();
					GRBLinExpr con6 = new GRBLinExpr();
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != i && m != i) {
								D6.set(cntr, getCol(x[i][k][m][j][r]), 1);
								con6_0.addTerm(1, x0[i][k][m][j][r]);
								con6.addTerm(1, x[i][k][m][j][r]);
							}
						}
					}
					D6.set(cntr, getLocCol(y[i]),M);
					d6.set(cntr, 0, M);
					con6_0.addTerm(M, y0[i]);
					con6.addTerm(M, y[i]);
					OP.addConstr(con6_0, GRB.LESS_EQUAL, M, "u6_" + i + "_" + j + "_" + r);
//					SP.addConstr(con6, GRB.LESS_EQUAL, M, "u6_" + i + "_" + j + "_" + r);
				}
			}
		}
		
		// Constraint 7
		Matrix D7 = new Matrix(N*(N-1)*(R+1)/2, nVar);
		Matrix d7 = new Matrix(N*(N-1)*(R+1)/2, 1, M);
		Matrix U7 = new Matrix(N*(N-1)*(R+1)/2, 1, 500);		
		cntr = 0; 		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= R; r++) {
					
					GRBLinExpr con7_0 = new GRBLinExpr();
					GRBLinExpr con7 = new GRBLinExpr();
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != j && m != j) {
								D7.set(cntr, getCol(x[i][k][m][j][r]), 1);
								con7_0.addTerm(1, x0[i][k][m][j][r]);
								con7.addTerm(1, x[i][k][m][j][r]);
							}
						}
					}
					D7.set(cntr,getLocCol(y[j]),M);
					d7.set(cntr, 0, M);
					con7_0.addTerm(M, y0[j]);
					con7.addTerm(M, y[j]);
					OP.addConstr(con7_0, GRB.LESS_EQUAL, M, "u7_" + i + "_" + j + "_" + r);
//					SP.addConstr(con7, GRB.LESS_EQUAL, M, "u7_" + i + "_" + j + "_" + r);
				}
			}
		}
		
		// Constraint 8
		Matrix D8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		Matrix U8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 500);
		cntr = 0; 
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= Math.pow(2, D) - 2; r++) {

					GRBLinExpr con8_0 = new GRBLinExpr();
					GRBLinExpr con8 = new GRBLinExpr();
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != i && k != j) {
								D8.set(cntr, getCol(x[i][k][m][j][r]), 1);
								if (r==0 && i == 0 && j == 1) System.out.println(x[i][k][m][j][r].get(GRB.StringAttr.VarName) + ": CoEf="+1+" - row="+cntr+" - col="+ getCol(x[i][k][m][j][r]));
								con8_0.addTerm(1, x0[i][k][m][j][r]);
								con8.addTerm(1, x[i][k][m][j][r]);
							}
						}
					}
					
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							D8.set(cntr, getCol(x[i][k][m][j][2 * r + 1]), -1);
							if (r==0 && i == 0 && j == 1) System.out.println(x[i][k][m][j][2 * r + 1].get(GRB.StringAttr.VarName) + ": CoEf= -1"+" - row="+cntr+" - col="+ getCol(x[i][k][m][j][2 * r +1]));
							con8_0.addTerm(-1, x0[i][k][m][j][2 * r + 1]);  // Left child node
							con8.addTerm(-1, x[i][k][m][j][2 * r + 1]);  // Left child node
						}
					}
					OP.addConstr(con8_0, GRB.LESS_EQUAL, 0, "u8_" + i + "_" + j + "_" + r);
//					SP.addConstr(con8, GRB.LESS_EQUAL, 0, "u8_" + i + "_" + j + "_" + r);
					cntr++;
				}
			}
		}
		
		// Constraint 9
		Matrix D9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		Matrix U9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 500);
		cntr = 0; 		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= Math.pow(2, D) - 2; r++) { 

					GRBLinExpr con9_0 = new GRBLinExpr();
					GRBLinExpr con9 = new GRBLinExpr();
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (m!=i && m!=j){
								D9.set(cntr, getCol(x[i][k][m][j][r]),1);
								con9_0.addTerm(1, x0[i][k][m][j][r]);
								con9.addTerm(1, x[i][k][m][j][r]);
							}
						}
					}
					
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {							
							D9.set(cntr, getCol(x[i][k][m][j][2 * r + 2]), -1);
							con9_0.addTerm(-1, x0[i][k][m][j][2 * r + 2]);
							con9.addTerm(-1, x[i][k][m][j][2 * r + 2]);
						}
					}
					OP.addConstr(con9_0, GRB.LESS_EQUAL, 0, "u9_" + i + "_" + j + "_" + r);
//					SP.addConstr(con9, GRB.LESS_EQUAL, 0, "u9_" + i + "_" + j + "_" + r);
					cntr++;
				}
			}
		}
		
		// Constraint 10
		Matrix D10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, M);
		Matrix U10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 500);
		cntr = 0; 		
		for (int i=0;i<N;i++){
			for (int j=i+1;j<N;j++){
				for (int k=0;k<N;k++){
					for (int r=0;r<=Math.pow(2, D) - 2;r++){
						
						GRBLinExpr con10_0 = new GRBLinExpr();
						GRBLinExpr con10 = new GRBLinExpr();
						for (int s:BinaryTree.getLeftChildren(r, D)){
							for (int m=0; m<N; m++){
								D10.set(cntr, getCol(x[i][k][m][j][s]), 1);
								D10.set(cntr, getCol(x[i][m][k][j][s]), 1);
								con10_0.addTerm(1, x0[i][k][m][j][s]);
								con10_0.addTerm(1, x0[i][m][k][j][s]);
								con10.addTerm(1, x[i][k][m][j][s]);
								con10.addTerm(1, x[i][m][k][j][s]);
							}
						}
						for (int m=0; m<N; m++){
							D10.set(cntr, getCol(x[i][k][m][j][r]), M);
							con10_0.addTerm(M , x0[i][k][m][j][r]);
							con10.addTerm(M , x[i][k][m][j][r]);
						}
						OP.addConstr(con10_0, GRB.LESS_EQUAL, M, "u10_" + i + "_" + j + "_" + k + "_" + r);
//						SP.addConstr(con10, GRB.LESS_EQUAL, M, "u10_" + i + "_" + j + "_" + k + "_" + r);
					}
				}
			}
		}
		
		// Constraint 11
		Matrix D11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, M);
		Matrix U11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 500);
		cntr = 0; 		
		for (int i=0;i<N;i++){
			for (int j=i+1;j<N;j++){
				for (int m=0;m<N;m++){
					for (int r=0;r<=Math.pow(2, D) - 2;r++){

						GRBLinExpr con11_0 = new GRBLinExpr();
						GRBLinExpr con11 = new GRBLinExpr();
						for (int s : BinaryTree.getRightChildren(r, D)) {
							for (int k = 0; k < N; k++) {
								D11.set(cntr, getCol(x[i][k][m][j][s]), 1);
								D11.set(cntr, getCol(x[i][m][k][j][s]), 1);
								con11_0.addTerm(1, x0[i][k][m][j][s]);
								con11_0.addTerm(1, x0[i][m][k][j][s]);
								con11.addTerm(1, x[i][k][m][j][s]);
								con11.addTerm(1, x[i][m][k][j][s]);
							}
						}
						for (int k = 0; k < N; k++) {
							D11.set(cntr, getCol(x[i][k][m][j][r]), M);
							con11_0.addTerm(M, x0[i][k][m][j][r]);
							con11.addTerm(M, x[i][k][m][j][r]);
						}
						OP.addConstr(con11_0, GRB.LESS_EQUAL, M, "u11_" + i + "_" + j + "_" + m + "_" + r);
//						SP.addConstr(con11, GRB.LESS_EQUAL, M, "u11_" + i + "_" + j + "_" + m + "_" + r);
					}
				}
			}
		}
		/*
		 * Lagrangian Relaxation Algorithm starts:
		 */
		// Obtaining Upper Bound
		double startTime = System.currentTimeMillis();
		SP.optimize();
//		printSol(SP);
		LB = obtainLB(OP, SP);
		UB = SP.get(GRB.DoubleAttr.ObjVal);
		
		ArrayList<Matrix> ds_List = new ArrayList<Matrix>();
		ArrayList<Matrix> Ds_List = new ArrayList<Matrix>();
		ArrayList<Matrix> Us_List = new ArrayList<Matrix>();
		ds_List.add(d6);ds_List.add(d7);ds_List.add(d8);ds_List.add(d9);ds_List.add(d10);ds_List.add(d11);
		Ds_List.add(D6);Ds_List.add(D7);Ds_List.add(D8);Ds_List.add(D9);Ds_List.add(D10);Ds_List.add(D11);	
		Us_List.add(U6);Us_List.add(U7);Us_List.add(U8);Us_List.add(U9);Us_List.add(U10);Us_List.add(U11);		
		ds = concatH(ds_List);
		Ds = concatH(Ds_List);
		Us = concatH(Us_List);
		System.out.println("miu: " + miu + " - Itr" + k + ": LB=" + LB.value + " - UB= " + UB);
		cntr = 0;
		double currentGap = GRB.INFINITY;
		while(currentGap > gap && (k<10 || LB.value!=bestLB)){
			cntr = dissectEpsilon(cntr, 15);
//			miu = updateMiu(miu, k, N);
			miu = updateMiuC(Ds, ds);		// Update Miu
			Us = updateU(SP, miu, Us, ds, Ds);		// Update Lagrangian multipliers
			updateObjCoeffs(SP, Us, ds, Ds);	// Update obj fun coefficients
			SP.optimize();
			UB = SP.get(GRB.DoubleAttr.ObjVal);
			cntr = updateUB(cntr, UB);
			LB = obtainLB(OP, SP);
			updateLB(LB.value);   // update best LB
			k++;	// Counter update
			currentGap = Math.abs((UB - LB.value)/LB.value);
			System.out.println("miu: " + miu + " - Itr" + k + ": LB=" + LB.value + " - UB= " + UB + " - Gap = " + currentGap + " - sol: " + printSol2(SP));
		}	
		List<String> selectedLocations = new ArrayList<String>();
		for (int i = 1 ; i <= N ; i++){
			GRBVar var = SP.getVar(SP.get(GRB.IntAttr.NumVars)-i);
			if (var.get(GRB.DoubleAttr.X) > 0)
				selectedLocations.add(var.get(GRB.StringAttr.VarName));					
		}
		System.out.println("Lower Bound: " + LB.value + "  -  Upper Bound: " + UB);
		
		/*
		 *  Snyder Daskin Variable Fixing procedure
		 */
		/*BenefitList benefits = computeBenefits(SP, x, y);
		fixVar2(benefits, y, UB, LB.value);
		System.out.println(fixTo0);
		System.out.println(fixTo1);
		for (int i = 0 ; i < N ; i++){
			if (y[i].get(GRB.DoubleAttr.X)==1){
				SP.get(GRB.DoubleAttr.ObjVal) + benefits.
			}else{
				
			}
		}*/
		
		/*
		 *  Variable fixing procedure
		 */
		// Getting the incumbent value
		/*setBinaries(OP);
		OP.optimize();
		for (int i = 0 ; i < N ; i++)
			fixVar(y0[i], y[i].get(GRB.DoubleAttr.X), "fix"+i);
		OP.optimize();
		printSol(OP);
		z = OP.get(GRB.DoubleAttr.ObjVal);
		z = LB.value;
		// Putting variables to be fixed into a list
		setN1N0(SP.getVars());
		varsTobeFixed(N1, N0, z);
		System.out.println("Done");
		
		// Removing fixing constraints
		for (int i = N ; i > 0 ; i--){
			int k = OP.get(GRB.IntAttr.NumConstrs) - 1;
			OP.remove(OP.getConstr(k));
			OP.update();
		}*/
		
		/*
		 * B&B procedure 2 
		 */
		BBNode rootNode = new BBNode(0, N, Us, LB, UB, selectedLocations);
		bestNode = rootNode;
		List<GRBVar> unfixedVars;
		unexploredNodes.add(rootNode);
		unfixedVars = getUnfixedVars(y, rootNode);
		rootNode.varToFix = getBranchVar(x, unfixedVars);
		System.out.println("y" + rootNode.varToFix.get(GRB.StringAttr.VarName));
		
		while(!unexploredNodes.isEmpty()){
			
			BBNode parent = unexploredNodes.get(0);
			System.out.println("Node index: " + parent.ind);
			System.out.print("Fixed to 1: ");
			for (GRBVar var: parent.varsFixedTo1){
				System.out.print(var.get(GRB.StringAttr.VarName));
			}
			System.out.println();
			System.out.print("Fixed to 0: ");
			for (GRBVar var: parent.varsFixedTo0){
				System.out.print(var.get(GRB.StringAttr.VarName));
			}
			System.out.println();
			
			if (!parent.fathom){
				
				BBNode rightChild = new BBNode(2*parent.ind+2, N, parent.Us, Ds, ds, OP, SP, bestLB, parent.varsFixedTo1, parent.varsFixedTo0, parent.varToFix, true);
				rightChild.updateFathom(P, N);
				if (!rightChild.fathom){
					unfixedVars = getUnfixedVars(y, rightChild);
					rightChild.varToFix = getBranchVar(x, unfixedVars);
//					updateLB(rightChild);				
					unexploredNodes.add(0, rightChild);
				}else{
					BBNodeList.add(rightChild);
				}
				LR_Main.relaxVar(SP, rightChild.noOfFixedVars);
				printNode(rightChild);
				
				
				BBNode leftChild = new BBNode(2*parent.ind+1, N, parent.Us, Ds, ds, OP, SP, bestLB, parent.varsFixedTo1, parent.varsFixedTo0, parent.varToFix, false);
				leftChild.updateFathom(P, N);
				if (!leftChild.fathom){
					unfixedVars = getUnfixedVars(y, leftChild);
					leftChild.varToFix = getBranchVar(x, unfixedVars);
//					updateLB(rightChild);
					unexploredNodes.add(0, leftChild);
				}else{
					BBNodeList.add(leftChild);
				}				
				LR_Main.relaxVar(SP, leftChild.noOfFixedVars);
				printNode(leftChild);
			}
			updateBestNode(parent);
//			exploredNodes.add(parent);
			BBNodeList.add(parent);
			unexploredNodes.remove(parent);
		}
		System.out.println(BBNodeList.peek().selectedLocations);
		System.out.println("Elapsed Time: " + (System.currentTimeMillis() - startTime));
		System.out.println("B&B done!");
		
		
		// Adding rows
		/*while(true){
			Matrix rhs = d8.minus(D8.times(X.transpose()));
			for (int i=0 ; i<rhs.getRowDimension() ; i++){
				if (rhs.get(i,0)<0){
					expr = new GRBLinExpr();
					System.out.println();
					for (int j = 0 ; j<nVar ; j++){
						expr.addTerm(D8.get(i, j),OP.getVar(varHash.));
						if (D8.get(i, j) != 0) System.out.println(D8.get(i, j) + "," + OP.getVar(j).get(GRB.StringAttr.VarName));

					}
					OP.addConstr(expr, GRB.LESS_EQUAL, 0, null);
				}
			}
			OP.optimize();
			printSol(OP);
			X = getVarsMatrix(OP);
		}*/
	}

}
