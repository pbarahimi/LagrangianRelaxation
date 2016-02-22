/*
 * Lower Bound Class
 */


import java.util.ArrayList;
import java.util.List;


public class LB {
	public double value;	// the value of the objective function to the original problem as Lower Bound
	public List<String> vars = new ArrayList<String>();	// names of the location variable with vlaue 1
	
	public LB(double value){
		this.value = value;
	}
	
	public LB(){
		
	}
}
