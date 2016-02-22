import gurobi.GRB;
import gurobi.GRBException;
import gurobi.GRBVar;

class VarIndices{
	protected int i;
	protected int j;
	protected int k;
	protected int m;
	protected int r;
	
	public VarIndices(GRBVar var) throws GRBException{
		String[] indices = var.get(GRB.StringAttr.VarName).split("_");
		this.i = Integer.parseInt(indices[0]);
		this.k = Integer.parseInt(indices[1]);
		this.m = Integer.parseInt(indices[2]);
		this.j = Integer.parseInt(indices[3]);
		this.r = Integer.parseInt(indices[4]);				
	}			
}