import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;

public class Buckets {
	int l_param;
	int tuples;
	int dims;
	int bucketSize; // = (tuples/l_param);
	Map<Short, Integer> sValues; //SA to frequences map.
	short[][] frequencies;	// [|SA|][2]
	short[][] map;			// [n][dims]
	short[][][] buckets;	// [l][(n/l)][dims]
	double cost;
	

	public Buckets(int l_param, int tuples, int dims, short[][] map,
				   short[][][] buckets) {
		this.l_param = l_param;
		this.tuples = tuples;
		this.dims = dims;
		this.map = map;
		this.buckets = buckets;
		bucketSize = (tuples/l_param);
	}

	public int getL() {
		return l_param;
	}

	public int getSize() {
		return tuples;
	}
	
	public void quickSortSA(int left, int right) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		int pivot = frequencies[ref][1];
		short[] temp2 = new short[2];
		while (i <= j) {
			//DESCENDING ORDER:
			while (pivot < (frequencies[i][1]))
				i++;
			while ((frequencies[j][1]) < pivot)
				j--;
			if (i <= j) {
				temp2=frequencies[i];
				frequencies[i]=frequencies[j];
				frequencies[j]=temp2;
				i++;
				j--;
			}
		}
		// recursion
		if (left < j)
			quickSortSA(left, j);
		if (i < right) {
			quickSortSA(i, right);
		}
	}
	
	public void quickSort(int left, int right, int SA) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		short pivot = map[ref][SA];
		int pfreq = sValues.get(pivot);
		short[] temp2 = new short[dims];
		while (i <= j) {
			//DESCENDING ORDER:
			while ( (pfreq < sValues.get(map[i][SA])) ||
				   ((pfreq == sValues.get(map[i][SA])) &&
				   (pivot < (map[i][SA]))) )
				i++;
			while ( (sValues.get(map[j][SA]) < pfreq) ||
				   ((sValues.get(map[j][SA]) == pfreq) &&
					((map[j][SA]) < pivot)) )
				j--;
			if (i <= j) {
				temp2=map[i];
				map[i]=map[j];
				map[j]=temp2;
				i++;
				j--;
			}
		}
		// recursion
		if (left < j)
			quickSort(left, j, SA);
		if (i < right)
			quickSort(i, right, SA);
	}
	
	public int getEmptyBucket(int[] empty){
		int best=0;
		for(int i=1; i<l_param; i++){
			if(empty[i] > empty[best])
				best = i;
		}
		return best;
	}
	
	public void printBuckets(){
		// Print tuples in buckets:
		for(int b=0; b<l_param; b++){
			System.out.println("bucket="+b);
			for(int t=0; t<bucketSize; t++){
				System.out.print("[");
				for(int a=0; a<dims; a++){
					System.out.print(" "+buckets[b][t][a]);
				}
				System.out.println(" ]");
			}
		}
	}
	
	public short[][][] bucketization(int SA){
		System.out.println("1. check if SA is in dims range");
		if ((SA<0)||(SA>=dims)) {
			//SA is out of range
			System.err.println("Error: "+SA+" is out of range");
			//SA = dims;
			System.exit(-1);
		}

		//count frequencies:
		System.out.println("2. count frequencies");
		sValues = new HashMap<Short, Integer>();
		int maxFreq = 1;
		for (int i = 0; i < tuples;i++) {
			if (sValues.containsKey(map[i][SA])) {
				sValues.put((map[i][SA]), sValues.get(map[i][SA]) + 1);
				if (maxFreq < (sValues.get(map[i][SA]))){
					maxFreq = (sValues.get(map[i][SA]));
				}
			}else{
				sValues.put((map[i][SA]), 1);
			}
		}
		
		//sort tuples wrt SA frequencies:
		System.out.println("3. sort tuples wrt SA frequencies");
		quickSort(0, tuples-1, SA);
		
		int frSize = sValues.size();
		/*
		frequencies = new short[frSize][2];
		int i=0;
		for (Map.Entry<Short, Integer> entry : sValues.entrySet()) {
			frequencies[i][0] = entry.getKey();
			frequencies[i][1] = entry.getValue();
		}
		sValues = null; //hopefully the garbage collector will free this.
		//System.gc(); //merely a suggestion...
		quickSortSA(0, frSize-1);
		*/
		
		//l-eligibility:
		System.out.println("4. l-eligibility");
		if (frSize < l_param){
			System.err.println("Error: There are less than l SA values");
			System.exit(-1);
		}
		if (maxFreq > bucketSize){ //(tuples / l_param)
			System.err.println("Error: max SA frequency("+maxFreq
							   +") is > than n/l ("+bucketSize+")");
			System.exit(-1);
		}
		
		//buckets:
		int offset = 0;
		int[] empty = new int[l_param];
		int freq = 0;
		int bucketOffset = 0;
		
		//first l most frequent SA values:
		System.out.println("5. first l most frequent SA values");
		for (int i=0; i<l_param; i++){
			freq = sValues.get(map[offset][SA]);
			/*
			System.out.println("SA="+(map[offset][SA])+"\tb="+i
							   +" fits="+freq+"/"+freq);
			*/
			for (int j=0; j<freq; j++){
				buckets[i][j] = map[offset + j];
			}
			empty[i] = bucketSize - freq;
			offset += freq;
		}
		
		//remaining SA values:
		System.out.println("6. remaining SA values");
		int fits = 0;
		int remains = 0;
		while(offset < tuples){
			freq = sValues.get(map[offset][SA]);
			remains = freq;
			while(remains > 0){
				int i = getEmptyBucket(empty);
				if(empty[i] < remains){
					fits = empty[i];
				}else{
					fits = remains;
				}
				bucketOffset = bucketSize - empty[i];
				for (int j=0; j<fits; j++){
					buckets[i][bucketOffset + j] = map[offset + j];
				}
				empty[i] -= fits;
				remains -= fits;
				offset += fits;
				/*
				System.out.println("SA="+(map[offset][SA])+"\tb="+i
								   +" fits="+fits+"/"+freq);
				*/
			}
		}
		
		//printBuckets(); //For debugging.
		
		System.out.println("7. return buckets");
		return buckets;
	}
}
