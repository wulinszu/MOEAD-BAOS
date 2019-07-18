
package moead_BAOS;

import jmetal.core.*;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.*;


public class MOEAD_IRA extends Algorithm {

	private int populationSize_;
	/**
	 * Stores the population
	 */
	private SolutionSet population_;
	/**
	 * Stores the values of the individuals
	 */
	private Solution[] savedValues_;

	private double[] utility_;
	private int[] frequency_;
	private double[] p;
	private int[] nc;
	
	double[] pn;

	/**
	 * Z vector (ideal point)
	 */
	double[] z_;
	
	double[] nz_;
	/**
	 * Lambda vectors
	 */

	double[][] lambda_;
	/**
	 * T: neighbour size
	 */
	int T_;
	/**
	 * Neighborhood
	 */
	int[][] neighborhood_;
	
	/**
	 * delta: probability that parent solutions are selected from neighbourhood
	 */
	double delta_;
	/**
	 * nr: maximal number of solutions replaced by each child solution
	 */
	int nr_;
	Solution[] indArray_;
	String functionType_;
	int evaluations_;
	int maxEvaluations;
	/**
	 * Operators
	 */
	Operator crossover_;
	Operator mutation_;

	String dataDirectory_;
	
	int runtime,fun;

	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public MOEAD_IRA(Problem problem) {
		super(problem);

		functionType_ = "_TCHE2";

	}
	
	public MOEAD_IRA(Problem problem, int runtime, int fun) {
		super(problem);

		functionType_ = "_TCHE2";
		this.runtime=runtime;
		this.fun = fun;

	} 

	public SolutionSet execute() throws JMException, ClassNotFoundException {
		

		evaluations_ = 0;
		maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations"))
				.intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize"))
				.intValue();
		dataDirectory_ = this.getInputParameter("dataDirectory").toString();
		
		T_ = 20;
		delta_ = 0.8;
		// nr_ = 2;

		population_ = new SolutionSet(populationSize_);
		savedValues_ = new Solution[populationSize_];
		utility_ = new double[populationSize_];
		frequency_ = new int[populationSize_];
		p = new double[populationSize_];
		pn = new double[T_];
		nc = new int[populationSize_];
		for (int i = 0; i < utility_.length; i++) {
			utility_[i] = 0;
			frequency_[i] = 0;
			p[i] = 0.5;
			nc[i] = 0;
		}
		indArray_ = new Solution[problem_.getNumberOfObjectives()];

		neighborhood_ = new int[populationSize_][T_];

		z_ = new double[problem_.getNumberOfObjectives()];
		nz_ 		  = new double[problem_.getNumberOfObjectives()];
		lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

		crossover_ = operators_.get("crossover"); // default: DE crossover
		mutation_ = operators_.get("mutation"); // default: polynomial mutation

		// STEP 1. Initialization
		// STEP 1.1. Compute euclidean distances between weight vectors and find
		// T
		initUniformWeight();
		initNeighborhood();

		// STEP 1.2. Initialize population
		initPopulation();

		// STEP 1.3. Initialize z_
		initIdealPoint();
		initNadirPoint();
		
		comp_p_neighbor();
		
		int gen = 0;
		// STEP 2. Update
		do {
			
			for(int n = 0; n < populationSize_; n++){
				
				double rand = PseudoRandom.randDouble();
				if (rand < p[n]) {
					
					frequency_[n]++;

					int type;
					double rnd = PseudoRandom.randDouble();

					// STEP 2.1. Mating selection based on probability
					if (rnd < delta_) 
					{
						type = 1; // neighborhood
					} else {
						type = 2; // whole population
					}
					Vector<Integer> p = new Vector<Integer>();
					
					//matingSelection(p, n, 2, type);
					
					matingSelection1(p, n, 5, type);

					// STEP 2.2. Reproduction
					Solution child;
					Solution[] parents = new Solution[5];
					parents[0] = population_.get(p.get(0));
					parents[1] = population_.get(p.get(1));
					parents[2] = population_.get(p.get(2));
					parents[3] = population_.get(p.get(3));
					parents[4] = population_.get(p.get(4));
					String DE_VARIANT=BAOS.execute(population_, populationSize_, problem_.getNumberOfObjectives(), n, neighborhood_, T_);
					crossover_.setParameter("DE_VARIANT", DE_VARIANT);

					// Apply DE crossover
					child = (Solution) crossover_.execute(new Object[] {
							population_.get(n), parents });

					// Apply mutation
					mutation_.execute(child);

					// Evaluation
					problem_.evaluate(child);

					evaluations_++;

					// STEP 2.3. Repair. Not necessary

					// STEP 2.4. Update z_
					updateReference(child);

					// STEP 2.5. Update of solutions
					updateProblem_original(child, n, type);
					
				}
			} // for

			gen++;
			if (gen % 20 == 0) {
				comp_pnc();
			}

		} while (evaluations_ < maxEvaluations);

		return population_;
	}
	
	public void print(String path,int[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");
	        bw.newLine();        
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  }
	
	public void comp_p_neighbor() {
		for(int i = 0; i < T_; i++){
			
			pn[i] = 0.05 + 0.95 * (1-1/(1+0.05*Math.exp(-20*((double)i/(double)20-0.7))));
		}
	}

	/**
	 * initUniformWeight
	 */
	public void initUniformWeight() {
		String dataFileName;
		dataFileName = "W" + problem_.getNumberOfObjectives() + "D_"
				+ populationSize_ + ".dat";

		// System.out.println(dataDirectory_);
		// System.out.println(dataDirectory_ + "/" + dataFileName);

		try {
			// Open the file
			FileInputStream fis = new FileInputStream(dataDirectory_ + "/"
					+ dataFileName);
			// FileInputStream fis = new FileInputStream(dataDirectory_);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);

			int numberOfObjectives = 0;
			int i = 0;
			int j = 0;
			String aux = br.readLine();
			while (aux != null) {
				StringTokenizer st = new StringTokenizer(aux);
				j = 0;
				numberOfObjectives = st.countTokens();
				while (st.hasMoreTokens()) {
					double value = (new Double(st.nextToken())).doubleValue();
					lambda_[i][j] = value;
					// System.out.println("lambda["+i+","+j+"] = " + value)
					// ;
					j++;
				}
				aux = br.readLine();
				i++;
			}
			br.close();
		} catch (Exception e) {
			System.err
					.println("initUniformWeight: failed when reading for file: "
							+ dataDirectory_ + "/" + dataFileName);
			e.printStackTrace();
		}

		// System.exit(0) ;
	} // initUniformWeight
	

	public void comp_pnc() {

		double f1, f2, maxUtility = 0;
		double e = 1.0 * java.lang.Math.pow(10, -50);
		double[][] distance = new double[populationSize_][populationSize_];
		int max_nc = 0;

		for (int i = 0; i < populationSize_; i++) {
			nc[i] = 0;
		}

		for (int n = 0; n < populationSize_; n++) {

			int minIndex = 0;
			for (int m = 0; m < populationSize_; m++) {// subproblem
				distance[n][m] = calculateDistance(population_.get(n),
						lambda_[m]);
				if (distance[n][m] < distance[n][minIndex]) {
					minIndex = m;
				}
			}
			nc[minIndex] = nc[minIndex] + 1;

			f1 = fitnessFunction(population_.get(n), lambda_[n]);
			f2 = fitnessFunction(savedValues_[n], lambda_[n]);
			double temp = (f2 - f1) / f2;
			utility_[n] = temp > 0 ? temp : 0;

			if (utility_[n] > maxUtility) {
				maxUtility = utility_[n];
			}
			savedValues_[n] = new Solution(population_.get(n));
		}

		for (int k = 0; k < populationSize_; k++) {
			if (nc[k] > max_nc) {
				max_nc = nc[k];
			}
		}

		for (int i = 0; i < populationSize_; i++) {
			p[i] = 0.98 * ((utility_[i] + e) / (maxUtility + e)) + 0.02
					* (1 - (double) nc[i] / (double) max_nc);
		}
	}
	
	public double calculateDistance(Solution individual, double[] lambda) {
		
		double scale;
		double distance;

		double[] vecInd  = new double[problem_.getNumberOfObjectives()];
		double[] vecProj = new double[problem_.getNumberOfObjectives()];
		
		// vecInd has been normalized to the range [0,1]
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			vecInd[i] = (individual.getObjective(i) - z_[i]) / (nz_[i] - z_[i]);

		scale = innerproduct(vecInd, lambda) / innerproduct(lambda, lambda);
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			vecProj[i] = vecInd[i] - scale * lambda[i];

		distance = norm_vector(vecProj);

		return distance;
		
	}

public double innerproduct(double[] vec1, double[] vec2) {
	double sum = 0;
	
	for (int i = 0; i < vec1.length; i++)
		sum += vec1[i] * vec2[i];
	
	return sum;
}

	public double norm_vector(double[] z) {
		double sum = 0;

		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			sum += z[i] * z[i];

		return Math.sqrt(sum);
	}

	/**
   *
   */
	public void initNeighborhood() {
		double[] x = new double[populationSize_];
		int[] idx = new int[populationSize_];

		for (int i = 0; i < populationSize_; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < populationSize_; j++) {
				x[j] = Utils.distVector(lambda_[i], lambda_[j]);
				// x[j] = dist_vector(population[i].namda,population[j].namda);
				idx[j] = j;
				// System.out.println("x["+j+"]: "+x[j]+
				// ". idx["+j+"]: "+idx[j]) ;
			} // for

			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, populationSize_, T_);
			// minfastsort(x,idx,population.size(),niche);

			System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
		} // for
	} // initNeighborhood

	/**
   *
   */
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			evaluations_++;
			population_.add(newSolution);
			savedValues_[i] = new Solution(newSolution);
		} // for
	} // initPopulation

	/**
   *
   */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			z_[i] = 1.0e+30;
			indArray_[i] = new Solution(problem_);
			problem_.evaluate(indArray_[i]);
			evaluations_++;
		} // for

		for (int i = 0; i < populationSize_; i++) {
			updateReference(population_.get(i));
		} // for
	} // initIdealPoint
	
	void initNadirPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			nz_[i] = -1.0e+30;

		for (int i = 0; i < populationSize_; i++)
			updateNadirPoint(population_.get(i));
	} // initNadirPoint

	/**
   *
   */
	public void matingSelection(Vector<Integer> list, int cid, int size,
			int type) {
		// list : the set of the indexes of selected mating parents
		// cid : the id of current subproblem
		// size : the number of selected mating parents
		// type : 1 - neighborhood; otherwise - whole population
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
				// p = population[cid].table[r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			// if (flag) list.push_back(p);
			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection
	
	public void matingSelection1(Vector<Integer> list, int cid, int size,
			int type) {
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		
		while (list.size() < size) {
			boolean flag = true;
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
				
				double rand = PseudoRandom.randDouble();
				if(rand > pn[r]){
					flag = false;
					//continue;
				}
				/*else{
					//System.out.println(rand);
					System.out.println(r);
					//System.out.println(pn[r]);
				}*/
				
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}
			if (flag) {
				list.addElement(p);
				//System.out.println(r);
			}
		}
	} // matingSelection
	

	public double calculatedis(Solution individual,double[] z){
		double dis=0;
		for(int i = 0;i<problem_.getNumberOfObjectives();i++){
			dis += (individual.getObjective(i)-z[i])*(individual.getObjective(i)-z[i]);
		}
		return Math.sqrt(dis);
	}

	public List<Integer> tour_selection(int depth) {

		// selection based on utility
		List<Integer> selected = new ArrayList<Integer>();
		List<Integer> candidate = new ArrayList<Integer>();

		for (int k = 0; k < problem_.getNumberOfObjectives(); k++)
			selected.add(k); // WARNING! HERE YOU HAVE TO USE THE WEIGHT
								// PROVIDED BY QINGFU (NOT SORTED!!!!)

		for (int n = problem_.getNumberOfObjectives(); n < populationSize_; n++)
			candidate.add(n); // set of unselected weights

		while (selected.size() < (int) (populationSize_ / 5.0)) {
			// int best_idd = (int) (rnd_uni(&rnd_uni_init)*candidate.size()),
			// i2;
			int best_idd = (int) (PseudoRandom.randDouble() * candidate.size());
			// System.out.println(best_idd);
			int i2;
			int best_sub = candidate.get(best_idd);
			int s2;
			for (int i = 1; i < depth; i++) {
				i2 = (int) (PseudoRandom.randDouble() * candidate.size());
				s2 = candidate.get(i2);
				// System.out.println("Candidate: "+i2);
				if (utility_[s2] > utility_[best_sub]) {
					best_idd = i2;
					best_sub = s2;
				}
			}
			selected.add(best_sub);
			candidate.remove(best_idd);
		}
		return selected;
	}

	/**
	 * 
	 * @param individual
	 */
	void updateReference(Solution individual) {
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			if (individual.getObjective(n) < z_[n]) {
				z_[n] = individual.getObjective(n);

				indArray_[n] = individual;
			}
		}
	} // updateReference
	
	void updateNadirPoint(Solution individual) {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (individual.getObjective(i) > nz_[i])
				nz_[i] = individual.getObjective(i);
		}
	} // updateNadirPoint

	/**
	 * @param individual
	 * @param id
	 * @param type
	 */
	void updateProblem(Solution indiv, int id, int type) {
		// indiv: child solution
		// id: the id of current subproblem
		// type: update solutions in - neighborhood (1) or whole population
		// (otherwise)
		int size;
		int time;

		time = 0;

		if (type == 1) {
			size = neighborhood_[id].length;
		} else {
			size = population_.size();
		}
		int[] perm = new int[size];

		Utils.randomPermutation(perm, size);

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1) {
				k = neighborhood_[id][perm[i]];
			} else {
				k = perm[i]; // calculate the values of objective function
								// regarding the current subproblem
			}
			double f1, f2;

			f1 = fitnessFunction(population_.get(k), lambda_[k]);
			f2 = fitnessFunction(indiv, lambda_[k]);

			if (f2 < f1) {
				population_.replace(k, new Solution(indiv));
				// population[k].indiv = indiv;
				time++;
			}
			// the maximal number of solutions updated is not allowed to exceed
			// 'limit'
			if (time >= nr_) {
				return;
			}
		}
	} // updateProblem

	
	void updateProblem_original(Solution indiv, int id, int type) {
		
		double maxu = 0;
		int k = 0;
		for (int i = 0; i < populationSize_; i++) {
			double f1, f2, uti;
			f1 = fitnessFunction(population_.get(i), lambda_[i]);
			f2 = fitnessFunction(indiv, lambda_[i]);
			
			uti = (f1 - f2) / f1;
			
			if (uti > maxu) {
				maxu = uti;
				k = i;
			}
		}
		
		if(maxu > 0){
			population_.replace(k, new Solution(indiv));
		}
	}
	
	
	/**
	 * @param individual
	 * @param lambda
	 * @return
	 */

	double fitnessFunction(Solution individual, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				double diff = Math.abs(individual.getObjective(i) - z_[i]);

				double feval;
				if (lambda[i] == 0) {
					feval = 0.000001 * diff;
				} else {
					feval = diff * lambda[i];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			fitness = maxFun;
		} else if (functionType_.equals("_TCHE2")) {
			double maxFun = -1.0e+30;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				double diff = Math.abs(individual.getObjective(i) - z_[i]);

				double feval;
				if (lambda[i] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[i];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			fitness = maxFun;
		} else if (functionType_.equals("_PBI")) {
			double theta; // penalty parameter
			theta = 2.0;

			// normalize the weight vector (line segment)
			double nd = norm_vector(lambda);
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
				lambda[i] = lambda[i] / nd;

			double[] realA = new double[problem_.getNumberOfObjectives()];
			double[] realB = new double[problem_.getNumberOfObjectives()];

			// difference between current point and reference point
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realA[n] = (individual.getObjective(n) - z_[n]);

			// distance along the line segment
			double d1 = Math.abs(innerproduct(realA, lambda));

			// distance to the line segment
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realB[n] = (individual.getObjective(n) - (z_[n] + d1
						* lambda[n]));
			double d2 = norm_vector(realB);

			fitness = d1 + theta * d2;
		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type "
					+ functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation

}
