package moead_BAOS;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.LZ09.LZ09_F1;
import jmetal.problems.LZ09.LZ09_F2;
import jmetal.problems.LZ09.LZ09_F3;
import jmetal.problems.LZ09.LZ09_F4;
import jmetal.problems.LZ09.LZ09_F5;
import jmetal.problems.LZ09.LZ09_F6;
import jmetal.problems.LZ09.LZ09_F7;
import jmetal.problems.LZ09.LZ09_F8;
import jmetal.problems.LZ09.LZ09_F9;
import jmetal.problems.M2M.MOP1;
import jmetal.problems.M2M.MOP2;
import jmetal.problems.M2M.MOP3;
import jmetal.problems.M2M.MOP4;
import jmetal.problems.M2M.MOP5;
import jmetal.problems.M2M.MOP6;
import jmetal.problems.M2M.MOP7;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.CEC2009_UF1;
import jmetal.problems.cec2009Competition.CEC2009_UF10;
import jmetal.problems.cec2009Competition.CEC2009_UF2;
import jmetal.problems.cec2009Competition.CEC2009_UF3;
import jmetal.problems.cec2009Competition.CEC2009_UF4;
import jmetal.problems.cec2009Competition.CEC2009_UF5;
import jmetal.problems.cec2009Competition.CEC2009_UF6;
import jmetal.problems.cec2009Competition.CEC2009_UF7;
import jmetal.problems.cec2009Competition.CEC2009_UF8;
import jmetal.problems.cec2009Competition.CEC2009_UF9;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.fastHypervolume.wfg.Front;
import jmetal.qualityIndicator.fastHypervolume.wfg.WFGHV;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.comparators.DominanceComparator;

public class MOEAD_BAOS_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	static String savePF = "Y";
	static String savePS = "Y";
	static String saveMETRICS = "Y";
	static String saveTable_IGD = "Y";
	static String saveTable_HV = "Y";
	static String PFFILE;
	static String PSFILE;
	static String METRICSFILE;
	
	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;//java文件输出流，创建文件流
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;//OutputStreamWriter是字符流通向字节流的桥梁 
	      BufferedWriter bw      = new BufferedWriter(osw)        ;//缓冲区               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");//写到缓冲区
	        bw.newLine(); //换行       
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	
	public static void main(String[] args) throws JMException, SecurityException, IOException, ClassNotFoundException {
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("MOEAD.log");
		logger_.addHandler(fileHandler_);
		String AlgorithmName="MOEAD_DE_S4";

		for (int fun = 10; fun <= 28; fun++) {
			int runtimes = 51;
			double[] IGDarray = new double[runtimes];
			double[] Hypervolume = new double[runtimes];
			// double Execuion_time = 0.0;

			Problem problem = null; // The problem to solve
			Algorithm algorithm; // The algorithm to use
			Operator crossover; // Crossover operator
			Operator mutation; // Mutation operator

			QualityIndicator indicators = null; // Object to get quality
												// indicators
			HashMap<String, Double> parameters; // Operator parameters
			
			/*	create files	*/
			String  DATAFILE 	  = "H:/2019_7_5_BAOS_test/"+AlgorithmName;	createFolder(DATAFILE);
			if (savePF.equals("Y"))
			{
				PFFILE 	  = DATAFILE + "/"+ "PF";	createFolder(PFFILE);
				//PFFILE 	  = PFFILE + "/"+ algName;	createFolder(PFFILE);
			}
			if (savePS.equals("Y"))
			{
				PSFILE 	  = DATAFILE + "/"+"PS";	createFolder(PSFILE);
				//PSFILE 	  = PSFILE + "/"+ algName;	createFolder(PSFILE);
			}
			if (saveMETRICS.equals("Y"))
			{
				METRICSFILE 	  = DATAFILE + "/"+"METRICS";	createFolder(METRICSFILE);
				//METRICSFILE 	  = METRICSFILE + "/"+ algName;	createFolder(METRICSFILE);
				//String METRICSPATH = METRICSFILE + "/" + problem.getName() + "M" + problem.getNumberOfObjectives();	
				//indicators.printVectorToFile(METRICSPATH, new String[]{"IGD"}, false);
			}

			System.out.println("============================================================================");
			double IGD_sum = 0;
			double HV_sum = 0;
			for (int run = 0; run < runtimes; run++) {

				if (fun == 1) {
					// problem = new WFG1("Real", 4, 20, 2);
					problem = new WFG1("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG1_605.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG1_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 2) {
					// problem = new WFG2("Real", 4, 20, 2);
					problem = new WFG2("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG2_111.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG2_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 3) {
					// problem = new WFG3("Real", 4, 20, 2);
					problem = new WFG3("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG3_right.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG3_2D.txt");
				}
				if (fun == 4) {
					// problem = new WFG4("Real", 4, 20, 2);
					problem = new WFG4("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG4_1181.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG4_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 5) {
					// problem = new WFG5("Real", 4, 20, 2);
					problem = new WFG5("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG5_694.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG5_2D.txt");
				}
				if (fun == 6) {
					// problem = new WFG6("Real", 4, 20, 2);
					problem = new WFG6("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\wfg6_166.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG6_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 7) {
					// problem = new WFG7("Real", 4, 20, 2);
					problem = new WFG7("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG7_2435.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG7_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 8) {
					// problem = new WFG8("Real", 4, 20, 2);
					problem = new WFG8("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG8_201.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG8_2D.txt");
				}
				if (fun == 9) {
					// problem = new WFG9("Real", 4, 20, 2);
					problem = new WFG9("Real", 8, 20, 2);
					// indicators = new QualityIndicator(problem,
					// "E:\\new_multiobjective\\jMetal\\Pareto_front\\WFG9_2591.txt");
					indicators = new QualityIndicator(problem,
							"H:\\2019_4_11_experiments\\truePF\\WFG2D\\WFG9_2D.txt");
				} // problem = new WFG1("Real");
				if (fun == 10) {
					problem = new CEC2009_UF1("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF1.dat");
				}
				if (fun == 11) {
					problem = new CEC2009_UF2("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF2.dat");
				}
				if (fun == 12) {
					problem = new CEC2009_UF3("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF3.dat");
				}
				if (fun == 13) {
					problem = new CEC2009_UF4("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF4.dat");
				}
				if (fun == 14) {
					problem = new CEC2009_UF5("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF5.dat");
				}
				if (fun == 15) {
					problem = new CEC2009_UF6("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF6.dat");
				}
				if (fun == 16) {
					problem = new CEC2009_UF7("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF7.dat");
				}
				if (fun == 17) {
					problem = new CEC2009_UF8("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF8.dat");
				}
				if (fun == 18) {
					problem = new CEC2009_UF9("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF9.dat");
				}
				if (fun == 19) {
					problem = new CEC2009_UF10("Real");
					indicators = new QualityIndicator(problem, "H:\\truePareto_front\\UF\\UF10.dat");
				}
				// -----------------------------
				if (fun == 20) {
					problem = new LZ09_F1("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F1.txt");
				}
				if (fun == 21) {
					problem = new LZ09_F2("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F2.txt");
				}
				if (fun == 22) {
					problem = new LZ09_F3("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F3.txt");
				}
				if (fun == 23) {
					problem = new LZ09_F4("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F4.txt");				
					}
				if (fun == 24) {
					problem = new LZ09_F5("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F5.txt");			
				}
				if (fun == 25) {
					problem = new LZ09_F6("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F6.txt");		
				}
				if (fun == 26) {
					problem = new LZ09_F7("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F7.txt");	
				}
				if (fun == 27) {
					problem = new LZ09_F8("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F8.txt");	
				}
				if (fun == 28) {
					problem = new LZ09_F9("Real");
					indicators = new QualityIndicator(problem,
							"H:\\truePareto_front\\LZ09_F\\LZ09_F9.txt");				
				}
				if (fun == 29) {
					problem = new MOP1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP1.dat");
				}
				if (fun == 30) {
					problem = new MOP2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP2.dat");
				}
				if (fun == 31) {
					problem = new MOP3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP3.dat");
				}
				if (fun == 32) {
					problem = new MOP4("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP4.dat");
				}
				if (fun == 33) {
					problem = new MOP5("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP5.dat");
				}
				if (fun == 34) {
					problem = new MOP6("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP6.dat");
				}
				if (fun == 35) {
					problem = new MOP7("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\pf_MOP7.dat");
				}
				if (fun == 36) {
					problem = new ZDT1("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT1_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 37) {
					problem = new ZDT2("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT2_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 38) {
					problem = new ZDT3("Real");

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT3_269.txt");
				}
				if (fun == 39) {
					problem = new ZDT4("Real", 10);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT4_501.txt");
				} // problem = new WFG1("Real");
				if (fun == 40) {
					problem = new ZDT6("Real", 10);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\ZDT6_774.txt");
				}
				if (fun == 41) {
					problem = new DTLZ1("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ1.pf");
				}
				if (fun == 42) {
					problem = new DTLZ2("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ2.pf");
				} // problem = new WFG1("Real");
				if (fun == 43) {
					problem = new DTLZ3("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ3.pf");
				} // problem = new WFG1("Real");
				if (fun == 44) {
					problem = new DTLZ4("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ4.pf");
				}
				if (fun == 45) {
					problem = new DTLZ5("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ5.txt");
				} // problem = new WFG1("Real");
				if (fun == 46) {
					problem = new DTLZ6("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ6.txt");
				} // problem = new WFG1("Real");
				if (fun == 47) {
					problem = new DTLZ7("Real", 10, 3);

					indicators = new QualityIndicator(problem,
							"E:\\new_multiobjective\\jMetal\\Pareto_front\\DTLZ7.pf");
				}

				algorithm = new MOEAD_DE(problem);

				Solution referencePoint = new Solution(problem);
				int popSize;
				int evaSize;
				if (fun >= 1 && fun <= 9) {// WFG1-9
					popSize = 100;
					evaSize = 25000;

					referencePoint.setObjective(0, 2.2);
					referencePoint.setObjective(1, 4.4);
				} else if (fun >= 10 && fun <= 16) {// UF1-7
					popSize = 600;
					evaSize = 300000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
				} else if (fun >= 17 && fun <= 19) {// UF8-10
					popSize = 1000;
					evaSize = 300000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
					referencePoint.setObjective(2, 1.1);
				} else if (fun >= 20 && fun <= 24) {// F1-F5
					popSize = 600;
					evaSize = 150000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
				} else if (fun == 25) {// F6
					popSize = 1000;
					evaSize = 300000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
					referencePoint.setObjective(2, 1.1);
				} else if (fun >= 26 && fun <= 28) {// F7-F9 (fun >= 26 && fun
													// <= 28)
					popSize = 600;
					evaSize = 150000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
				} else if (fun >= 29 && fun <= 33) {// MOP1-5
					popSize = 100;
					evaSize = 300000;

				} else if (fun >= 36 && fun <= 40) {// ZDT1-5
					popSize = 100;
					evaSize = 25000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
				} else if (fun >= 41 && fun <= 47) {// DTLZ1-7
					popSize = 500;
					evaSize = 100000;

					referencePoint.setObjective(0, 1.1);
					referencePoint.setObjective(1, 1.1);
				} else {
					popSize = 300;
					evaSize = 300000;
				}

				algorithm.setInputParameter("populationSize", popSize);
				algorithm.setInputParameter("maxEvaluations", evaSize);
				algorithm.setInputParameter("dataDirectory", "H:\\truePareto_front\\Weight\\");

				/*
				 * algorithm.setInputParameter("T", 20);
				 * algorithm.setInputParameter("delta", 0.9);
				 * //algorithm.setInputParameter("delta", 0.8);
				 * algorithm.setInputParameter("nr", 2);
				 */

				parameters = new HashMap<String, Double>();
				parameters.put("CR", 1.0);
				parameters.put("F", 0.5);
				parameters.put("K", 0.5);
				// crossover =
				// CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover_MAB",parameters);
				crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover_BAOS", parameters);
				parameters = new HashMap<String, Double>();
				parameters.put("probability", 1.0 / problem.getNumberOfVariables());
				parameters.put("distributionIndex", 20.0);
				mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
				algorithm.addOperator("crossover", crossover);
				algorithm.addOperator("mutation", mutation);

				// long initTime = System.currentTimeMillis();
				long initTime = System.currentTimeMillis();
				SolutionSet population = algorithm.execute(); // Execute the
																// algorithm
				IGDarray[run]=((Double)indicators.getCECIGD(population));
				SolutionSet s = new SolutionSet(population.size());
				Comparator c = new DominanceComparator();
				for(int j=0; j<population.size(); j++) {
					int v = c.compare(population.get(j), referencePoint);
					if(v == -1) {
						s.add(population.get(j));
					}
				}
				if(s.size() <= 2) {
					Hypervolume[run] = 0.0;
				}
				else {
					Front front = new Front(s.size(), referencePoint.getNumberOfObjectives(), s);					
					WFGHV hv = new WFGHV(referencePoint.getNumberOfObjectives(), s.size(), referencePoint);
					Hypervolume[run] = hv.getHV(front);
				}				
				
				long estimatedTime = System.currentTimeMillis() - initTime;
				
				//System.out.println(problem.getName() + "	M:"+ problem.getNumberOfObjectives() + "	Run:" + Integer.toString(run) 
				//	+ "	TimeUsed:" + estimatedTime + "ms" );
				//System.out.println(problem.getName() + "	M:"+ problem.getNumberOfObjectives() + "	Run:" + Integer.toString(run) 
				//	+ "	TimeUsed:" + estimatedTime + "ms" + "	IGD:" + ((Double)indicators.getIGD(population)).toString()+ 
				//	"	HV:" + ((Double)indicators.getHypervolume(population)).toString());
					
				//	"  HV:no");
				
				if (savePF.equals("Y"))
				{
					String PFPATH = PFFILE + "/" + problem.getName()+ "M" + problem.getNumberOfObjectives() + "R" + Integer.toString(run);
					population.printObjectivesToFile(PFPATH);
				}
				if (savePS.equals("Y"))
				{
					String PSPATH = PSFILE + "/" + problem.getName()+ "M" + problem.getNumberOfObjectives() + "R" + Integer.toString(run);
					population.printVariablesToFile(PSPATH);
				}
				
				if (saveMETRICS.equals("Y"))
				{
					String METRICSPATH = METRICSFILE + "/" + problem.getName() + "M" + problem.getNumberOfObjectives()+"_IGD.txt";					
					indicators.printVectorToFile(METRICSPATH, 
						new Double[]{IGDarray[run]}, true);
						//new Double[]{indicators.getIGD(population)}, true);
					String METRICSPATH2 = METRICSFILE + "/" + problem.getName() + "M" + problem.getNumberOfObjectives()+"_HV.txt";					
					indicators.printVectorToFile(METRICSPATH2, 
						new Double[]{Hypervolume[run]}, true);
						//new Double[]{indicators.getIGD(population)}, true);
							
				}
		
				//System.out.println("run_"+run+",IGD="+IGDarray[run]+",HV="+Hypervolume[run]);
				System.out.println("run_"+run+",IGD="+IGDarray[run]);
				
				IGD_sum=IGD_sum+IGDarray[run];
				HV_sum=HV_sum+Hypervolume[run];
				//System.out.println();
				//System.out.println();
			}// for 30 runtimes
			//将所有测试实例的平均IGD值存入一个文本中
			if (saveTable_IGD.equals("Y"))
			{
				String METRICSPATH =  "H:\\2019_7_5_BAOS_test\\"+AlgorithmName+ "\\Table_IGD.txt";	
				if(fun==1){
					indicators.printVectorToFile(METRICSPATH,new Double[]{IGD_sum/runtimes}, false);
				}else{
					indicators.printVectorToFile(METRICSPATH,new Double[]{IGD_sum/runtimes}, true);
				}
				
			}
			//将所有测试实例的平均HV值存入一个文本中
			if (saveTable_HV.equals("Y"))
			{
				String METRICSPATH = "H:\\2019_7_5_BAOS_test\\"+AlgorithmName + "\\Table_HV.txt";	
				if(fun==1){
					indicators.printVectorToFile(METRICSPATH,new Double[]{HV_sum/runtimes}, false);
				}else{
					indicators.printVectorToFile(METRICSPATH,new Double[]{HV_sum/runtimes}, true);
				}
				
			}
			
			System.out.println(problem.getName() + "	M:"+ problem.getNumberOfObjectives() + 
					"	the avergal IGD:"+IGD_sum/runtimes);
			System.out.println(problem.getName() + "	M:"+ problem.getNumberOfObjectives() + 
					"	the avergal HV:"+HV_sum/runtimes);	
		}// for finish a problem (30 times)

	} // main
	
	private static void createFolder(String str) {
		boolean success = (new File(str)).mkdir();

		/*if (success) {
			System.out.println("Directory: " + str + " created");
		} else {
			System.out.println("Directory: " + str + " NOT created");
		}*/
	}

} 
