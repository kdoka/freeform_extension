package generators;

//package parsers;
import java.util.Random;

/**
 * 
 * @author Katerina
 *
 */
public class Stats {
	
	double[] probs;
	double sum;

	public Stats(int size, double theta) {
		probs = getProbs(size, -theta);
	}

	// random number generator
	private Random rand = new Random();

	public int getRandInt (int N) {
		return rand.nextInt(N);
	}
	public double getRand () {
		return rand.nextDouble();
	}

	public double[] getProbs(int size, double theta){
		double[] nums = new double[size];
		for (int i = 0; i < nums.length; i++) {
			nums[i] = i + 1;
		}
		double[] probs = new double[nums.length];
		for (int i = 0; i < probs.length; i++) {
			if (nums[i] == 0) { probs[i] = 0; }
			else { 
				probs[i] = Math.pow(nums[i],theta); 
			}
		}
		return probs;
	}

	
	public int select (double[] probs) {
		// make array of probabilities

		// sum probabilities
		sum = 0;
		for (int i = 0; i < probs.length; i++) {
			sum += probs[i];
		}
		// obtain random number in range [0, sum]
		double r = sum * getRand();
		// subtract probs until result negative
		// no of iterations gives required index
		int i;
		for (i = 0; i < probs.length; i++) {
			r -= probs[i];
			if (r < 0) { break; }
		}
		return i;
	}
	
	public double accumProb (int stop){
		double accum = 0;
		for (int i = 0; i < stop; i++) {
			accum += probs[i];
		}
		return accum/sum;
	}

	// select item using Zipf's law
	// parameter is size of ranked array
	// returns index in [0, array size - 1]
	public int zipfGen () {
		// get index using special case of power law
		return select(this.probs);
	}

	
}