package ie.errity.pd.genetic;

import ie.errity.pd.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.BitSet;
import java.util.Random;
import java.util.ArrayList;

/**
 * Provides a means of evolving a population of 
 *{@link  ie.errity.pd.Prisoner Prisoners} 
 *via a genetic algorithm
 * @author	Andrew Errity 99086921
 * @author      Sherri Goings (modified 4/4/2017)
 */
public class Breeder extends JPanel
{
    private Prisoner curPopulation[];
    private double mutateP, crossP; //mutate, cross probabilities
    private int selection, selParam, seed; // selection method, any associated selection parameter, random number seed
    private int popSize;
    private Random rand;
	
    /**
     *Create a new Genetic Breeder
     */
    public Breeder()
    {	
	//defaults
	mutateP = .001;
	crossP = .95;
	selection = 0;
	selParam = 1;
	seed = -1;
	rand = new Random(); //uses current time as seed
    }
	
    /**
     *Create a new Genetic Breeder using {@link  ie.errity.pd.Rules Rules} given
     *@param r1 parameters for the breeder
     */
    public Breeder(Rules r1){
	mutateP = r1.getMutateP();
	crossP = r1.getCrossP();
	selection = r1.getSelection();
	selParam = r1.getSelectionParameter();
	seed = r1.getSeed();
	if (seed==-1)
	    rand = new Random();
	else
	    rand = new Random(seed); //uses seed set in rules for repeatability
    }

	
    /**
     *Breeds the next generation (panmictic mating) of an array of 
     *{@link  ie.errity.pd.Prisoner Prisoners} 
     *@param c	initial population (raw fitness of population must be calcualted previously)
     *@return the next generation
     */
    public Prisoner[] Breed(Prisoner[] c){
		curPopulation = c;	//population to breed
		popSize = curPopulation.length;
		Prisoner Selected[] = new Prisoner[popSize]; // parent pop after selection 
		
			
		// Select parents for next gen 
		// ***ADD CODE (make separate functions) to perform each of the types of selection***
		// selection is an int that determines the method to use
		// selParam is an int that can be used as a parameter for a selection method if required
		
		// One silly selection method which uses rand and the selParam is given below
		// Be sure to use "rand" class variable here as your random number generator (i.e. don't create a new one!)

		// selection method 0, use parameter as a threshhold percentage.  If a random number is less than param, 
		// select the best individual from the population.  Otherwise select a random individual.
		// In this method a param of 0 gives random selection, and a param of 100 gives best-wins-all selection.
		if (selection == 0) {
			// find index of most fit individual
			double maxFit = 0;
			int indexBest = 0;
			for (int i=0; i<popSize; i++) {
				if (curPopulation[i].getScore() > maxFit) {
					maxFit = curPopulation[i].getScore();
					indexBest = i;
				}
			}
			
			// select according to description above for this method
			for (int i=0; i<popSize; i++) {
				int selIndex = 0;
				if (rand.nextInt(100) < selParam){
					selIndex = indexBest;
				}
				else{
					selIndex = rand.nextInt(popSize);
				}
				Selected[i] = (Prisoner)curPopulation[selIndex].clone();
			}
		}
		
		// selection method 1 implements optional elitism, fitness proportional & stochastic universal sampling
		else if(selection == 1){
			sortByFitness(curPopulation);
			
			// Optinal elitism
			ArrayList<Prisoner> elites = null;
			if(selParam > 0){
				if(selParam > popSize){
					System.out.println("Error: selParam > popSize");
					System.exit(0);
				}else{
					elites = new ArrayList<Prisoner>();
					for(int i = 0; i < selParam; i++){
						elites.add(curPopulation[i]);
					}
				}
			}
			
			// Sigma scaling
			double[] scaledFitnesses = scaleFitnessValues(curPopulation);
			
			// Fitness proportional & stochastic universal sampling
			Prisoner[] newSelected = FPSUS(curPopulation, scaledFitnesses);
			
			int selIndex = 0;
			for(int i = 0; i < elites.size(); i++){
				Selected[selIndex++] = (Prisoner)elites.get(i).clone();
			}
			for(int i = selIndex; i < popSize; i++){
				Selected[i] = (Prisoner)newSelected[i].clone();
			}
		}
		
		// selection method 2 implements tournament selection
		else if(selection == 2){
			int selIndex = 0;
			double k = 0.75;
			while(selIndex < popSize){
				
				// generate two random indices that are not equal to each other
				int idx1 = 0;
				int idx2 = 0;
				while(idx1 == idx2){
					idx1 = getRandomNumber(0,popSize-1);
					idx2 = getRandomNumber(0,popSize-1);
				}
				
				Prisoner p1 = curPopulation[idx1];
				Prisoner p2 = curPopulation[idx2];
				Prisoner tournament_selected = null;
				double r = Math.random();
				
				// select high fitness
				if(r < k){
					tournament_selected = (p1.getScore() > p2.getScore()) ? p1 : p2;
					//System.out.println("Higher chosen: " + p1.getScore() + ", " + p2.getScore());
				}
				
				// select lower fitness
				else{
					tournament_selected = (p1.getScore() > p2.getScore()) ? p2 : p1;
					//System.out.println("Lower chosen: " + p1.getScore() + ", " + p2.getScore());
				}
				Selected[selIndex++] = (Prisoner)tournament_selected.clone();
			}
		}else {  // any other selection method fill pop with always cooperate
			for (int i=0; i<popSize; i++)
			Selected[i] = new Prisoner("ALLC");
		}

				
		//Crossover & Mutate each pair of selected parents	
		BitSet Offspring[] = new BitSet[2];  // temporarily holds 2 children during crossover/mutation
		for (int d=0; d<popSize; d+=2) {
			// in case of odd population, just mutate and replace last individual
			if (d+1 >= popSize) {
				Offspring[0] = Genetic.mutate(Selected[d].getStrat(), mutateP, rand);
				Selected[d] = new Prisoner(Offspring[0]);
			}
			else {			
				if(rand.nextDouble() <= crossP) //Cross Over
					Offspring = Genetic.crossover(Selected[d].getStrat(),Selected[d+1].getStrat(), rand);
				else{
					Offspring[0] = (BitSet)Selected[d].getStrat().clone();
					Offspring[1] = (BitSet)Selected[d+1].getStrat().clone();
				}
					
				//Mutation
				Offspring[0] = Genetic.mutate(Offspring[0],mutateP, rand);
				Offspring[1] = Genetic.mutate(Offspring[1],mutateP, rand);
					
				//Replacement - we are done with parents d & d+1, so just replace with children without
				// creating an entire new array
				Selected[d] = new Prisoner(Offspring[0]);
				Selected[d+1] = new Prisoner(Offspring[1]);
			}
		}
		// pass on children pop to be parents of next gen
		curPopulation = Selected;
		repaint();	//update display (if any)
		return curPopulation; //return the bred population
    }
	
	/**
	 * Returns a random number given a range
	 */
	private int getRandomNumber(int min, int max){
		Random rand = new Random();
		return rand.nextInt((max-min)+1) + min;
	}
	
	/**
     * Given a list of scaled fitness values,
     * select offsprings using fitness proportional & stochastic universal sampling
     */
    private Prisoner[] FPSUS(Prisoner[] curPopulation, double[] scaledFitness){
        int size = curPopulation.length;
		double total_fitness = 0;
		
        for(int i = 0; i < scaledFitness.length; i++){
            total_fitness += scaledFitness[i];
        }
        
        double ptr = total_fitness/size;
        double start = Math.random();
        double[] ptrs = new double[size];
        for(int i = 0; i < size; i++){
            ptrs[i] = start + i*ptr;
        }
        return RWS(curPopulation, scaledFitness, ptrs);
    }
	
	private Prisoner[] RWS(Prisoner[] population, double[] scaledFitness, double[] points){
        Prisoner[] selected = new Prisoner[points.length];
        double[] runningSum = new double[points.length];
        for(int i = 0; i < runningSum.length; i++){
            if(i == 0){
                runningSum[i] = scaledFitness[i];
            }else{
                runningSum[i] = runningSum[i-1] + scaledFitness[i];
            }
        }
        
        for(int i = 0; i < points.length; i++){
            double p = points[i];
            int index = 0;
            while((runningSum[index] < p)){
                index++;
                if (index == points.length-1){
                break;
                }
            }
            selected[i] = (Prisoner)population[index].clone();
        }
        return selected;
    }
	
	/**
	 * Given a list of prionsers, return a list of scaled fitness values
	 */
	private double[] scaleFitnessValues(Prisoner[] c){
		double[] scaledFitness = new double[c.length];
		double fitTotal = 0;
		for (int i = 0; i < c.length; i++){
			fitTotal += c[i].getScore();
		}
    
		double fitMean = fitTotal/c.length;
		double variance = 0;
		double stdDev = 0;
    
		for (int i = 0; i < c.length; i++){
			variance += Math.pow((c[i].getScore() - fitMean),2);		
		}
    
		variance = variance / c.length;
		stdDev = Math.sqrt(variance);
    
		for (int i = 0; i < c.length; i++){
			if(stdDev == 0){
				scaledFitness[i] = 1;
			}else{
				scaledFitness[i] = 1 + (c[i].getScore()-fitMean)/(2*stdDev);
			}
		}
		return scaledFitness;
    }
	
	/**
	 * Sort population by fitness
	 */
	private void sortByFitness(Prisoner[] c){
		quickSort(c, 0, c.length-1);
	}
	
	private void quickSort(Prisoner[] c, int low, int high){
		if(low < high){
			int p = partition(c, low, high);
			quickSort(c, low, p-1);
			quickSort(c, p+1, high);
		}
	}
	
	private int partition(Prisoner[] c, int low, int high){
		int pivot = c[high].getScore();
		int wall = low;
		for(int i = low; i < high; i++){
			if(c[i].getScore() > pivot){
				Prisoner temp = c[i];
				c[i] = c[wall];
				c[wall] = temp;
				wall++;
			}
		}
		Prisoner temp = c[wall];
		c[wall] = c[high];
		c[high] = temp;
		return wall;
	}
	
    /**
     *Responsible for updating the graphical representation
     */
    public void paintComponent(Graphics g) 
    {
        super.paintComponent(g); //paint background

	//Get limits of viewable area
      	Insets insets = getInsets();
      	int x0 = insets.left;
      	int y0 = insets.top;
	int currentWidth = getWidth() - insets.left - insets.right;
	int currentHeight = getHeight() - insets.top - insets.bottom;
	    
	//Display a series of rectangles, representing the players
	for(int i = 0; i < popSize; i++)
	    {
	    	g.setColor(curPopulation[i].getColor());	
	    	g.fillRect((x0*2)+((currentWidth/popSize)*(i)),(currentHeight/4)+y0,(currentWidth/popSize),currentHeight/2);
	    }
    }
	
    /**
     *Set the {@link  ie.errity.pd.Rules Rules}
     *@param r new {@link  ie.errity.pd.Rules Rules}
     */
    public void setRules(Rules r)
    {
	mutateP = r.getMutateP();
	crossP = r.getCrossP();
	selection = r.getSelection();
	selParam = r.getSelectionParameter();
	seed = r.getSeed();
	if (seed > -1)
	    rand.setSeed(seed);
    }
	
    /**
     *Reset the breeder
     */
    public void clear()
    {
	popSize = 0;
	repaint();
    }
}
	
