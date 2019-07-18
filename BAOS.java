package moead_BAOS;

import java.util.Comparator;
import java.util.Vector;

import MOEAD_Test2.Utils;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.util.PseudoRandom;
import jmetal.util.comparators.DominanceComparator;

public class BAOS {
	
	//population_: the whole population
	//populationSize_: the size of population_
	//m: the number of objectives
	//n: the index of current selected solution in population
	//T: the size of neighborhood of each solution
	//neighborhood: the neighborhood indexs of all solutions
	//parents: parent solutions for recombination
	
	public static String execute(SolutionSet population_, int populationSize_, int m, int n, int[][] neighborhood, int T) {
		String DE_VARIANT=null;
		//Pareto-assisted Criterion: label==1非支配解（注重exploitation）,label==2被支配（注重exploration）
		int label=1;
		Comparator dominance_ = new DominanceComparator();
		for (int k = 0; k < (populationSize_ - 1); k++) {
			int flag = dominance_.compare(population_.get(n), population_.get(k));
			if (flag==1){
				label=2;
				break;
			}
		}		

		// STEP 2.2. Reproduction
		//Crowding-assisted Criterion: strategy_flag==1注重exploitation,strategy_flag==2注重exploration
		int strategy_flag;
		switch(label){
		case 1:
			//Crowding-assisted Criterion
			strategy_flag=CrowdingCompetition(population_,populationSize_,neighborhood,n,m,T);
			if (strategy_flag==1){
				DE_VARIANT="DE/rand/1";
			}else{
				DE_VARIANT="DE/current-to-rand/1";
			}
			break;
		case 2:
			//Crowding-assisted Criterion
			strategy_flag=CrowdingCompetition(population_,populationSize_,neighborhood,n,m,T);
			if (strategy_flag==1){
				DE_VARIANT="DE/rand/2";
			}else{
				DE_VARIANT="DE/current-to-rand/2";
			}
			break;
		}//switch
		return DE_VARIANT;
	}
	
	public static int CrowdingCompetition(SolutionSet population_,int populationSize_, int [][] neighborhood_, int n,int m,int T_){
		int strategy_flag;
		SolutionSet SS=new SolutionSet(T_);
		SolutionSet SSS=new SolutionSet(T_);
		int[] permutation1 = new int[populationSize_];
		Utils.randomPermutation(permutation1, populationSize_);
		int temp=permutation1[0];
		for (int i=0;i<T_;i++){
			SS.add(population_.get(neighborhood_[n][i]));
			SSS.add(population_.get(neighborhood_[temp][i]));
		}
		
		double crowd1=density(population_.get(temp),SSS,m,T_);
		double crowd=density(population_.get(n),SS,m,T_);
		if (crowd>crowd1){
			strategy_flag=1;
		}
		else if (crowd<crowd1){
			strategy_flag=2;
		}
		else{
			double rnd1 = PseudoRandom.randDouble();
			if (rnd1<0.5){
				strategy_flag=1;
			}else{
				strategy_flag=2;
			}
		}
		return strategy_flag;
	}
	
	public static int CrowdingCompetition_n(SolutionSet population_,int populationSize_, int [][] neighborhood_, int n,int m,int T_){
		int strategy_flag;
		SolutionSet SS=new SolutionSet(T_);
		SolutionSet SSS=new SolutionSet(T_);
		int[] permutation1 = new int[T_];
		Utils.randomPermutation(permutation1, T_);
		int temp=neighborhood_[n][permutation1[0]];
		for (int i=0;i<T_;i++){
			SS.add(population_.get(neighborhood_[n][i]));
			SSS.add(population_.get(neighborhood_[temp][i]));
		}
		
		double crowd1=density(population_.get(temp),SSS,m,T_);
		double crowd=density(population_.get(n),SS,m,T_);
		if (crowd>crowd1){
			strategy_flag=1;
		}
		else if (crowd<crowd1){
			strategy_flag=2;
		}
		else{
			double rnd1 = PseudoRandom.randDouble();
			if (rnd1<0.5){
				strategy_flag=1;
			}else{
				strategy_flag=2;
			}
		}
		return strategy_flag;
	}
	
	public static double density(Solution individual,SolutionSet S,int m,int T){
		double dist=0.0;
		double sum;
		double [] d=new double[T];
		for (int i=0;i<T;i++){
			sum=0;
			for(int j=0;j<m;j++){
				sum+=Math.pow(individual.getObjective(j)-S.get(i).getObjective(j),2);
			}
			d[i]=Math.sqrt(sum);
		}
		
		for (int j=1;j<T;j++){
			if (d[j]==0){
				d[j]=1.0e-10;
			}
			dist=dist+1.0/d[j];
		}

		double crowd=(double)(T-1)/dist;
		return crowd;
	}
	
	public static String execute_random12(){
		String DE_VARIANT=null;
		double rnd = PseudoRandom.randDouble();
		if (rnd<0.25){
			DE_VARIANT="DE/rand/1";
		}
		else if (rnd>=0.25&&rnd<0.5){
			DE_VARIANT="DE/current-to-rand/1";
		}
		else if (rnd>=0.5&&rnd<0.75){
			DE_VARIANT="DE/rand/2";
		}
		else{
			DE_VARIANT="DE/current-to-rand/2";
		}
		return DE_VARIANT;
	}
	
	public static String execute_random1(SolutionSet population_, int populationSize_, int m, int n, int[][] neighborhood, int T){
		String DE_VARIANT=null;
		int strategy_flag;
		double rnd = PseudoRandom.randDouble();
		if (rnd<0.5){
			strategy_flag=CrowdingCompetition(population_,populationSize_,neighborhood,n,m,T);
			if (strategy_flag==1){
				DE_VARIANT="DE/rand/1";
			}else{
				DE_VARIANT="DE/current-to-rand/1";
			}
		}
		else{
			strategy_flag=CrowdingCompetition(population_,populationSize_,neighborhood,n,m,T);
			if (strategy_flag==1){
				DE_VARIANT="DE/rand/2";
			}else{
				DE_VARIANT="DE/current-to-rand/2";
			}
		}
		
		return DE_VARIANT;
	}
	
	public static String execute_random2(SolutionSet population_, int populationSize_,int n){
		String DE_VARIANT=null;
		//Pareto-assisted Criterion: label==1非支配解（注重exploitation）,label==2被支配（注重exploration）
		int label=1;
		Comparator dominance_ = new DominanceComparator();
		for (int k = 0; k < (populationSize_ - 1); k++) {
			int flag = dominance_.compare(population_.get(n), population_.get(k));
			if (flag==1){
				label=2;
				break;
			}
		}		

		if (label==1){
			double rnd = PseudoRandom.randDouble();
			if (rnd<0.5){
				DE_VARIANT="DE/rand/1";
			}
			else{
				DE_VARIANT="DE/current-to-rand/1";
			}
		}
		else{
			double rnd1 = PseudoRandom.randDouble();
			if (rnd1<0.5){
				DE_VARIANT="DE/rand/2";
			}
			else{
				DE_VARIANT="DE/current-to-rand/2";
			}
		}
		return DE_VARIANT;
	}

}
