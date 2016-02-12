import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import Jama.Matrix;
import gurobi.*;

public class LR_test {
	private static GRBModel model;
	private static Matrix C;
	private static HashMap<String, Integer> varHash = new HashMap<String, Integer>();
	private static int k = 0;
	private static double miu = 0.5;
	private static double epsilon = 1;
	private static double UB = 14;
	private static Matrix Us;
	private static Matrix Ds;
	private static Matrix ds;
	
	// Variable fixing attributes
	private static double z; // incumbent value for variable fixing
	private static List<GRBVar> N0;
	private static List<GRBVar> N1;
	private static List<GRBVar> fixTo1;
	private static List<GRBVar> fixTo0;
	
	/**
	 * prints the solution
	 * @param model
	 * @throws GRBException
	 */
	private static void printSol() throws GRBException{
		for (GRBVar var : model.getVars()) System.out.println(var.get(GRB.StringAttr.VarName) + " : " + var.get(GRB.DoubleAttr.X));
	}
	
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
	private static void updateObjCoeffs(Matrix U, Matrix d, Matrix D) throws GRBException{
		double[][] temp =  C.getArrayCopy();
		Matrix updatedC = new Matrix(temp);
		double constant = 0;
		updatedC = updatedC.minus(U.transpose().times(D));			
		constant += (U.transpose().times(d)).get(0, 0);
		for (GRBVar var: model.getVars()){
			int index = varHash.get(var.get(GRB.StringAttr.VarName));
			var.set(GRB.DoubleAttr.Obj, updatedC.getArray()[0][index]);
		}
		model.set(GRB.DoubleAttr.ObjCon, constant);	
	}

	/**
	 * updates the step-size based on the iteration number
	 * @param miu
	 * @param k
	 * @return
	 */
	private static double updateMiu(){
		int itr = (int) Math.floor(k/4);
		miu = miu* Math.pow(0.5,itr);
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
	 * Updates the miu based on the  rule (c)
	 * @param model
	 * @param varHash
	 * @param UB
	 * @param currentU
	 * @param epsilon
	 * @return
	 * @throws GRBException 
	 */
	private static double updateMiuC(Matrix D, Matrix d) throws GRBException{
		Matrix X = toMatrix(model.getVars(), varHash);
		Matrix denuminatorMat = D.times(X.transpose());		
		denuminatorMat = d.minus(denuminatorMat);
		double denuminator = Math.pow(denuminatorMat.normF(),2);	 
		double output = epsilon*(model.get(GRB.DoubleAttr.ObjVal) - UB ) / denuminator;
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
	private static Matrix updateU(Matrix U, Matrix d, Matrix D) throws GRBException{
		Matrix X1 = toMatrix(model.getVars(), varHash);
		Matrix output = D.times(X1.transpose());
		
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
	
	private static Matrix toMatrix(GRBVar[] vars, HashMap<String,Integer> varHash) throws GRBException{
		Matrix output = new Matrix(1, vars.length);
		for (GRBVar var:vars){
			output.set(0, varHash.get(var.get(GRB.StringAttr.VarName)), var.get(GRB.DoubleAttr.X));;			
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
			col = varHash.get(var.get(GRB.StringAttr.VarName)); // column index of the variable
			String s = var.get(GRB.StringAttr.VarName);
			Matrix temp = new Matrix(Us.getRowDimension(), Us.getColumnDimension());
			temp = Ds.getMatrix(0, Ds.getRowDimension()-1, col, col);
			temp = Us.arrayTimes(temp);
			double d = sumElements(temp);
			double k = C.get(0, col);
			double crt = k-d;
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
			int col = varHash.get(var.get(GRB.StringAttr.VarName));
			Matrix temp = Ds.getMatrix(0, rowSize, col, col);
			temp = Us.arrayTimes(temp);
			double tmep2 = sumElements(temp);
			double temp2 = C.get(0, col);
			temp2 = temp2 - tmep2;
			statement2 += C.get(0, col) - sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));			
		}
		
		// Defining variables to be fixed to 0
		for (GRBVar var : N1){
			int col = varHash.get(var.get(GRB.StringAttr.VarName));
			criteria = statement1 + statement2 + C.get(0, col) - sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));
			if (criteria >= z) fixTo0.add(var);
			System.out.println("cri: "+ criteria);
			}
		// Defining variables to be fixed to 1
		double statement2_2;
		for (GRBVar var : N0){
			int col = varHash.get(var.get(GRB.StringAttr.VarName));
			statement2_2 = statement2 - C.get(0, col) + sumElements(Us.arrayTimes(Ds.getMatrix(0, rowSize, col, col)));
			criteria = statement1 + statement2_2;
			System.out.println("cri: "+ criteria);
			if (criteria >= z) fixTo1.add(var);
		}		
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
	
	public static void main(String[] arg) throws GRBException{
		GRBEnv env = new GRBEnv("RpHND_LR.log");
		model = new GRBModel(env);
		model.getEnv().set(GRB.IntParam.OutputFlag, 0);
		
		// initializing vars
		GRBVar x1 = model.addVar(0, 1, 4, GRB.BINARY, "x1");
		GRBVar x2 = model.addVar(0, 1, 5, GRB.BINARY, "x2");
		GRBVar x3 = model.addVar(0, 1, 6, GRB.BINARY, "x3");
		GRBVar x4 = model.addVar(0, 1, 7, GRB.BINARY, "x4");
		GRBVar[] vars = {x1,x2,x3,x4};
		model.update();
		
		// creating variables hashMap
		int temp = 0;
		for (GRBVar var : model.getVars()){
			varHash.put(var.get(GRB.StringAttr.VarName),temp++);
		}
		
		// obj function declaration
		GRBLinExpr expr = new GRBLinExpr();
		double[] coEffs = {4,5,6,7};
		expr.addTerms(coEffs, vars);
		model.setObjective(expr, GRB.MAXIMIZE);
		double[][] coEffs0 = {{4,5,6,7}};
		C = new Matrix(coEffs0);		
		
		// constraint 1
		double[][] coEffs1 = {{2,2,3,4}};
		Matrix D1 = new Matrix(coEffs1);
		double[][] rhs1 = {{7}};
		Matrix d1 = new Matrix(rhs1);
		Matrix u1 = new Matrix(1, 1);
		
		//constraint 2
		double[][] coEffs2 = {{1,-1,1,-1}};
		Matrix D2 = new Matrix(coEffs2);
		double[][] rhs2 = {{0}};
		Matrix d2 = new Matrix(rhs2);
		Matrix u2 = new Matrix(1, 1);
		
		/*
		 * Lagrangian Relaxation Algorithm starts:
		 */
		double secondLastObjVal = -1*GRB.INFINITY;
		model.optimize();
		double lastObjVal = model.get(GRB.DoubleAttr.ObjVal);
		List<Matrix> Ds_temp = new ArrayList<Matrix>();
		List<Matrix> ds_temp = new ArrayList<Matrix>();
		Ds_temp.add(D1); Ds_temp.add(D2);
		ds_temp.add(d1); ds_temp.add(d2);
		Ds = concatH(Ds_temp);
		ds = concatH(ds_temp);
		Us = new Matrix(2,1);
		while(miu>0.00001){
			miu = updateMiu();
//			miu = updateMiuC(Ds, ds);
			Us = updateU(Us, ds, Ds);	// Update the multipliers
			updateObjCoeffs(Us, ds, Ds);
			model.optimize();
			k = k + 1;
			printSol();
			/*for (int i = 0 ;i<Us.getRowDimension();i++){
				Matrix temp2 = Ds.getMatrix(i, i, 0 , Ds.getColumnDimension()-1);
				temp2 = temp2.times(toMatrix(model.getVars(), varHash).transpose());
				temp2 = ds.getMatrix(i,i,0,0).minus(temp2);
				System.out.println("Constraint"+i+" -> "+sumElements(temp2));
			}*/
			secondLastObjVal = lastObjVal;
			lastObjVal = model.get(GRB.DoubleAttr.ObjVal);
			System.out.println("miu: " + miu + " - u1:" + Us.get(0,0) + " - u2:" + Us.get(1, 0) + " - Itr" + k + ": " + lastObjVal);
		}
		printSol();
		
		/*
		 *  Variable fixing procedure
		 */
		GRBModel test = new GRBModel(env);
		/*double[][] D_temp = {{-1,0,-1,0,0,-1},{0,-1,0,-1,-1,0},{-1,-1,-1,0,0,0},{0,0,-1,0,-1,0}};
		double[][] d_temp = {{-1},{-1},{-1},{-1}};
		double[][] c_temp = {{-6,-7,-11,-5,-8,-5}};*/
		double[][] D_temp = {{1,0,1,0,0,1},{0,1,0,1,1,0},{1,1,1,0,0,0},{0,0,1,0,1,0}};
		double[][] d_temp = {{1},{1},{1},{1}};
		double[][] c_temp = {{6,7,11,5,8,5}};
		double[][] x_temp = {{1,0,0,0,1,0}};
		Ds = new Matrix(D_temp);
		ds = new Matrix(d_temp);
		C = new Matrix(c_temp);
		Matrix X = new Matrix(x_temp);
		Us = new Matrix(4,1);
		
		Us.set(0,0,4);
		Us.set(1,0,4);
		Us.set(2,0,3);
		Us.set(3,0,3);
		z = 14; // Incumbent value
		setN1N0(model.getVars());		
		varsTobeFixed(N1, N0, z);
		System.out.println("Fixed to 0:");
		for (GRBVar var : fixTo0) System.out.println(var.get(GRB.StringAttr.VarName));
		System.out.println("Fixed to 1:");
		for (GRBVar var : fixTo1) System.out.println(var.get(GRB.StringAttr.VarName));
		
		System.out.println("Done");
	}

}
