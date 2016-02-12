import java.util.HashMap;
import java.util.Map.Entry;

import gurobi.*;
import Jama.*;

public class LR_Main2 {
	// Lagrangian Relaxation algorithm attributes
	private static double gap = 0.001; // Termination criteria: the gap between the latest objVal and the second last. 
	private static double UB = -3430; // Upper Bound for the original problem.
	private static double epsilon = 2; // the constant used in the rule (c) of updating lagrangian multipliers update steps.
		
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
	 */
	
	/**
	 * Updates the objective function coefficients based on the Lagrangian Multipliers put in.
	 * @param model
	 * @param uHash
	 * @param U
	 * @param C
	 * @param d
	 * @param D
	 * @throws GRBException
	 */
	private static void updateObjCoeffs(GRBModel model, HashMap<String,Integer> varHash, Matrix U, Matrix C, Matrix d, Matrix D) throws GRBException{
		Matrix updatedC = C.minus(U.transpose().times(D));
		for (GRBVar var: model.getVars()){
			int index = varHash.get(var.get(GRB.StringAttr.VarName));
			var.set(GRB.DoubleAttr.Obj, updatedC.getArray()[0][index]);
		}
		double constant = (U.transpose().times(d)).get(0, 0);
		model.set(GRB.DoubleAttr.ObjCon, constant);
	}
	
	
	/**
	 * Updates the miu, Lagrangian multipliers and the objective function coefficients. 
	 * @param model
	 * @param uHash
	 * @param U
	 * @param C
	 * @param d
	 * @param D
	 * @throws GRBException
	 */
	private static void updateModel(GRBModel model, HashMap<String,Integer> varHash, Matrix U, Matrix C, Matrix d, Matrix D) throws GRBException{
		double miu = updateMiuC(model, varHash, UB, U, D, d, epsilon); // Update miu
		U = updateU(U, d, D, miu, model, varHash);	// Update the multipliers
		updateObjCoeffs(model, varHash, U, C, d, D);
	}
	
	/**
	 * Takes an array of Gurobi variables and return a 1xn matrix.
	 * @param vars
	 * @return 1xn Matrix
	 * @throws GRBException
	 */
	public static Matrix makeMat(GRBVar[] vars) throws GRBException{
		double[][] y = new double[0][vars.length];
		for (int i = 0; i<vars.length;i++){
			y[0][i] = vars[i].get(GRB.DoubleAttr.X);
		}
		return new Matrix(y);
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
	private static Matrix updateU(Matrix U, Matrix d, Matrix D, double miu, GRBModel model, HashMap<String,Integer> varHash) throws GRBException{
		Matrix X = toMatrix(model.getVars(), varHash);
		/*for (GRBVar var:model.getVars()){
			if (var.get(GRB.DoubleAttr.X)>0) System.out.println(var.get(GRB.StringAttr.VarName) + "-" + var.get(GRB.DoubleAttr.X));
		}*/
		Matrix U2 = D.times(X.transpose());
		U2.arrayTimesEquals(U);
		U2 = d.minus(U2);
		U2.timesEquals(miu);
		U2 = U.minus(U2);
		/*double[][] u_k = U2.getArray();
		for (int i=0;i<u_k.length;i++){
			for (int j=0;j<u_k[0].length;j++){
				u_k[i][j] = (u_k[i][j] > 0) ? u_k[i][j]:0;
			}
		}*/
		return U2;//new Matrix(u_k);		
	}
	
	private static Matrix toMatrix(GRBVar[] vars, HashMap<String,Integer> varHash) throws GRBException{
		Matrix output = new Matrix(1, vars.length);
		for (GRBVar var:vars){
			output.set(0, varHash.get(var.get(GRB.StringAttr.VarName)), var.get(GRB.DoubleAttr.X));;			
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
	
	/**
	 * Updates the miu based on the  rule (b)
	 * @param miu
	 * @param ro
	 * @return miu as double
	 */
	private static double updateMiuB(double miu, double ro){
		return miu * ro;
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
	private static double updateMiuC(GRBModel model, HashMap<String,Integer> varHash, double UB, Matrix U, Matrix D, Matrix d, double epsilon) throws GRBException{
		Matrix X = toMatrix(model.getVars(), varHash);
		Matrix denuminatorMat = d.minus((D.times(X.transpose())).arrayTimes(U));
		double denuminator = denuminatorMat.norm2();
//		System.out.println("denuminator: " + denuminator);		 
		double output = epsilon*(model.get(GRB.DoubleAttr.ObjVal) - UB ) / denuminator;
		return output;
	}
	
	public static void main(String[] args) throws GRBException {
		
		//Filling in the flows matrix asymmetrically 
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}
				
		// Building model
		GRBEnv env = new GRBEnv("RpHND_LR.log");
		GRBModel model = new GRBModel(env);
		model.getEnv().set(GRB.IntParam.OutputFlag, 0);
		
		// Creating variables
		GRBVar[][][][][] x = new GRBVar[N][N][N][N][R + 1];
		GRBVar[] y = new GRBVar[N];

		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						String varName = "x" + i + "_" + k + "_" + m + "_"
								+ j + "_0";
						x[i][k][m][j][0] = model.addVar(0.0, GRB.INFINITY, 0.0,
								GRB.BINARY, varName);						
					}
				}
			}
		}		

		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						for (int r = 1; r <= R; r++){
							String varName = "x" + i + "_" + k + "_" + m
									+ "_" + j + "_" + r;
							x[i][k][m][j][r] = model.addVar(0.0, GRB.INFINITY, 0.0,
									GRB.BINARY, varName);
						}
					}
				}
			}
		}

		for (int i = 0; i < N; i++) {
			y[i] = model.addVar(0, GRB.INFINITY, 0, GRB.BINARY, "y" + i);
		}
		

		// Integrate new variables
		model.update();

		// Updating HashMap of the variables - varHash
		HashMap<String,Integer> varHash = new HashMap<String, Integer>();
		int cntr = 0;
		for (GRBVar var: model.getVars()){
			varHash.put(var.get(GRB.StringAttr.VarName), cntr++);			
		}
		
		// Set objective function
		Matrix C = new Matrix(1,(int) ((Math.pow(N, 3)*(N-1)*(R+1)/2)+N));
		GRBLinExpr expr = new GRBLinExpr();
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						double CoEf = -1 * flows[i][j] * Cikmj(i, k, m, j) * (1 - Q(i,k,m,j));
						expr.addTerm(CoEf, x[i][k][m][j][0]);
						C.set(0, varHash.get(x[i][k][m][j][0].get(GRB.StringAttr.VarName)), CoEf);
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
							expr.addTerm(CoEf, x[i][k][m][j][r]);
							C.set(0, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), CoEf);
						}
					}
				}
			}
		}

		model.setObjective(expr, GRB.MAXIMIZE);
		
		// Adding constraints
		// Constraint 2
		GRBLinExpr con2 = new GRBLinExpr();
		for (int i = 0; i < N; i++) {
			con2.addTerm(1, y[i]);
		}
		model.addConstr(con2, GRB.EQUAL, P, "u2");

		// Constraint 3
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {

				GRBLinExpr con3 = new GRBLinExpr();
				for (int k = 0; k < N; k++) {
					for (int m = 0; m < N; m++) {
						con3.addTerm(1, x[i][k][m][j][0]);
					}
				}
				model.addConstr(con3, GRB.EQUAL, 1, "u3_" + i + "_" + j);
			}
		}

		// Constraint 4
		for (int r = 0; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int k = 0; k < N; k++) {

						GRBLinExpr con4 = new GRBLinExpr();
						for (int m = 0; m < N; m++) {
							con4.addTerm(1, x[i][k][m][j][r]);
						}
						con4.addTerm(-1, y[k]);
						model.addConstr(con4, GRB.LESS_EQUAL, 0, "u4_" + i + "_" + j + "_" + k + "_" + r);
					}
				}
			}
		}

		// Constraint 5
		for (int r = 0; r <= R; r++) {
			for (int i = 0; i < N; i++) {
				for (int j = i + 1; j < N; j++) {
					for (int m = 0; m < N; m++) {

						GRBLinExpr con5 = new GRBLinExpr();
						for (int k = 0; k < N; k++) {
							con5.addTerm(1, x[i][k][m][j][r]);
						}
						con5.addTerm(-1, y[m]);
						model.addConstr(con5, GRB.LESS_EQUAL, 0, "u5_" + i + "_" + j + "_" + m + "_" + r);
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
		Matrix U6 = new Matrix(N*(N-1)*(R+1)/2, 1);		
		D6 = new Matrix(N*(N-1)*(R+1)/2,nVar);
		HashMap<String, Integer> con6Hash = new HashMap<String, Integer>();
		cntr = 0;
		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= R; r++) {

					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != i && m != i) {
								D6.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), 1);
							}
						}
					}
					D6.set(cntr, varHash.get(y[i].get(GRB.StringAttr.VarName)),M);
					d6.set(cntr, 0, M);
					con6Hash.put("u6_" + i + "_" + j + "_" + r, cntr);
				}
			}
		}
		
		// Constraint 7
		Matrix D7 = new Matrix(N*(N-1)*(R+1)/2, nVar);
		Matrix d7 = new Matrix(N*(N-1)*(R+1)/2, 1, M);
		Matrix U7 = new Matrix(N*(N-1)*(R+1)/2, 1, 0);		
		HashMap<String, Integer> con7Hash = new HashMap<String, Integer>();
		cntr = 0; 
		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= R; r++) {

					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != j && m != j) {
								D7.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), 1);
							}
						}
					}
					D7.set(cntr,varHash.get(y[j].get(GRB.StringAttr.VarName)),M);
					d7.set(cntr, 0, M);
					con7Hash.put("u7_" + i + "_" + j + "_" + r, cntr++);					
				}
			}
		}
		
		// Constraint 8
		Matrix D8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		Matrix U8 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 1);
		HashMap<String, Integer> con8Hash = new HashMap<String, Integer>();
		cntr = 0; 
		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= Math.pow(2, D) - 2; r++) {

					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (k != i && k != j) {
								D8.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), 1);
							}
						}
					}

					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							D8.set(cntr, varHash.get(x[i][k][m][j][2 * r + 1].get(GRB.StringAttr.VarName)), -1);
						}
					}
					con8Hash.put("u8_" + i + "_" + j + "_" + r, cntr++);
				}
			}
		}
		
		// Constraint 9
		Matrix D9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		Matrix U9 = new Matrix((int) (N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 1);
		HashMap<String, Integer> con9Hash = new HashMap<String, Integer>();
		cntr = 0; 
		
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				for (int r = 0; r <= Math.pow(2, D) - 2; r++) { 

					GRBLinExpr con9 = new GRBLinExpr();
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							if (m!=i && m!=j){
								con9.addTerm(1, x[i][k][m][j][r]);
								D9.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)),1);
							}
						}
					}
					
					for (int k = 0; k < N; k++) {
						for (int m = 0; m < N; m++) {
							con9.addTerm(-1, x[i][k][m][j][2 * r + 2]);
							D9.set(cntr, varHash.get(x[i][k][m][j][2 * r + 2].get(GRB.StringAttr.VarName)), -1);
						}
					}
					con9Hash.put("u9_" + i + "_" + j + "_" + r, cntr++);
				}
			}
		}
		
		// Constraint 10
		Matrix D10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, M);
		Matrix U10 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		HashMap<String, Integer> con10Hash = new HashMap<String, Integer>();
		cntr = 0; 
		
		for (int i=0;i<N;i++){
			for (int j=i+1;j<N;j++){
				for (int k=0;k<N;k++){
					for (int r=0;r<=Math.pow(2, D) - 2;r++){
						
						for (int s:BinaryTree.getLeftChildren(r, D)){
							for (int m=0; m<N; m++){
								D10.set(cntr, varHash.get(x[i][k][m][j][s].get(GRB.StringAttr.VarName)), 1);
								D10.set(cntr, varHash.get(x[i][m][k][j][s].get(GRB.StringAttr.VarName)), 1);
							}
						}
						for (int m=0; m<N; m++){
							D10.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), M);
						}
						con10Hash.put("u10_" + i + "_" + j + "_" + k + "_" + r, M);
					}
				}
			}
		}
		
		// Constraint 11
		Matrix D11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), nVar);
		Matrix d11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, M);
		Matrix U11 = new Matrix((int) (N*N*(N-1)*(Math.pow(2, D) - 1)/2), 1, 0);
		HashMap<String, Integer> con11Hash = new HashMap<String, Integer>();
		cntr = 0; 
		
		for (int i=0;i<N;i++){
			for (int j=i+1;j<N;j++){
				for (int m=0;m<N;m++){
					for (int r=0;r<=Math.pow(2, D) - 2;r++){

						GRBLinExpr con11 = new GRBLinExpr();
						for (int s : BinaryTree.getRightChildren(r, D)) {
							for (int k = 0; k < N; k++) {
								D11.set(cntr, varHash.get(x[i][k][m][j][s].get(GRB.StringAttr.VarName)), 1);
								D11.set(cntr, varHash.get(x[i][m][k][j][s].get(GRB.StringAttr.VarName)), 1);
							}
						}
						for (int k = 0; k < N; k++) {
							D11.set(cntr, varHash.get(x[i][k][m][j][r].get(GRB.StringAttr.VarName)), M);
						}
						model.addConstr(con11, GRB.LESS_EQUAL, M, "u11_" + i + "_" + j + "_" + m + "_" + r);
						con11Hash.put("u11_" + i + "_" + j + "_" + m + "_" + r, M);
					}
				}
			}
		}
		
		/*
		 * Lagrangian Relaxation Algorithm starts:
		 */
		double secondLastObjVal = GRB.INFINITY;
		//solveSubProblem(model, varHash, U6, C, d6, D6);	// Update the objective function coefficients and solves the updated model to get new X's
		model.optimize();
		double lastObjVal = model.get(GRB.DoubleAttr.ObjVal);
		int cntr2 = 0;
		System.out.println("Itr" + cntr2 + ": " + lastObjVal);
		while(run(lastObjVal,secondLastObjVal)){
			updateModel(model, varHash, U6, C, d6, D6);
//			updateModel(model, varHash, U7, C, d7, D7);
//			updateModel(model, varHash, U8, C, d8, D8);
//			updateModel(model, varHash, U9, C, d9, D9);
//			updateModel(model, varHash, U10, C, d10, D10);
			model.optimize();
			secondLastObjVal = lastObjVal;
			lastObjVal = model.get(GRB.DoubleAttr.ObjVal);
			System.out.println("Itr" + cntr2 + ": " + lastObjVal);
		}
	}
}