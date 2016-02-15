import gurobi.GRB;
import gurobi.GRBException;
import gurobi.GRBVar;

import java.util.Collections;
import java.util.PriorityQueue;

public class BenefitList {
	
	PriorityQueue<Benefit> selectedLocs = new PriorityQueue<BenefitList.Benefit>();
	PriorityQueue<Benefit> unselectedLocs = new PriorityQueue<BenefitList.Benefit>(Collections.reverseOrder());
	
	public void add(GRBVar var, int index, double benefit) throws GRBException{
		if (var.get(GRB.DoubleAttr.X) == 1 )
			selectedLocs.add(new Benefit(index, benefit));
		else 
			unselectedLocs.add(new Benefit(index, benefit));
	}
	
	protected class Benefit implements Comparable<Benefit>{
		public int index;
		public double benefit;
		
		public Benefit(int index, double benefit){
			this.benefit = benefit;
			this.index = index;
		}

		@Override
		public int compareTo(Benefit other){
			double temp = this.benefit - other.benefit;
			if (temp>0) return 1;
			else if (temp<0) return -1;
			else return 0;
		}
		
		@Override
		public String toString(){
			return this.index + ":" + (int) this.benefit;
		}
	}
	

}
