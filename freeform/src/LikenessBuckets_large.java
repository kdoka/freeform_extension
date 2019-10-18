import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.ArrayList;
import java.io.IOException;
import java.io.FileWriter;

public class LikenessBuckets_large {
	double b_param;
	int tuples;
	int dims;
	int bucketSize; // c = (tuples/|B|);
	int buckNum; // |B| = Math.ceil(tuples/c);
	static Map<Short, ArrayList<Integer>> sTuples; //SA to tuples'index in map[].
	Map<Integer, ArrayList<Integer>> changeSA; //points of new SA per bucket.
	static short[] sArray;			// [|SA|]
	short[][] map;			// [n][dims]
	short[][][] buckets;	// [(n/c)][c][dims]
	double cost;
	static byte[] cardinalities = {79, 2, 17, 6, 9, 10, 83, 51};
	double threshold = 0.0;
	String inFile;
	

	public LikenessBuckets_large(double b_param, int tuples, int dims, short[][] map,
				   short[][][] buckets, double th, String inputFile) {
		this.b_param = (double)b_param;
		this.tuples = tuples;
		this.dims = dims;
		this.map = map;
		this.buckets = buckets;
		//bucketSize = 1;
		//buckNum=tuples;
		this.threshold = th;
		this.inFile = inputFile;//.substring(inputFile.lastIndexOf('/')+1, inputFile.lastIndexOf('.'));
	}

	public double getB() {
		return b_param;
	}

	public int getSize() {
		return tuples;
	}
	
	public static int GCD(int a, int b) {
		if (a==0 || b==0) {
			return Math.abs(a+b);
		}
		return GCD(b, a % b);
	}
	
	public static int gcd_all() {
		int g = 0;
		int i = sTuples.get(sArray[0]).size();
		int j = sTuples.get(sArray[1]).size();
		for(int ii=2; ii< sArray.length; ii++ ){
			g = GCD(i, j);
			i = g;
			j = sTuples.get(sArray[ii]).size();
		}
		return g;
	}
	
	public void quickSortSA_d(int left, int right) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		int pivot = sTuples.get(sArray[ref]).size();
		short temp2;
		while (i <= j) {
			//DESCENDING ORDER:
			while (pivot < sTuples.get(sArray[i]).size())
				i++;
			while (sTuples.get(sArray[j]).size() < pivot) 
				j--;
			if (i <= j) {
				temp2=sArray[i];
				sArray[i]=sArray[j];
				sArray[j]=temp2;
				i++;
				j--;
			}
		}
		// recursion
		if (left < j)
			quickSortSA_d(left, j);
		if (i < right) {
			quickSortSA_d(i, right);
		}
	}
	
	public void quickSortSA(int left, int right) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		int pivot = sTuples.get(sArray[ref]).size();
		short temp2;
		while (i <= j) {
			//ASCENDING ORDER:
			while (pivot > sTuples.get(sArray[i]).size())
				i++;
			while (sTuples.get(sArray[j]).size() > pivot)
				j--;
			if (i <= j) {
				temp2=sArray[i];
				sArray[i]=sArray[j];
				sArray[j]=temp2;
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
	
	//returns the most empty bucket
	public int getEmptyBucket(int[] empty){
		int best=0;
		for(int i=1; i<buckNum; i++){
			if(empty[i] > empty[best])
				best = i;
		}
		return best;
	}
	
	//returns the next most empty bucket, after currentB
	//or -1 if no other bucket is empty.
	public int getNextEmptyBucket(int[] empty, int currentB, int remains){
		int best=0;
				
		for(int i=1; i<buckNum; i++) {
			if (empty[i] > empty[best]){
				if ( (empty[i] < empty[currentB]) ||
					( (empty[i] == empty[currentB]) && (i>currentB) ) ) {
					best = i;
				}
			}
		}
		
		if ((empty[best] == 0) || (empty[best] > empty[currentB])) {
			//System.out.println(" - ");
			return -1;
		}
		
		if (best == currentB){
			//System.out.println(" -- ");
			return -1;
		}
		
		if ( (empty[currentB] >= remains) && (empty[best] < remains) ){ //*0.9
			//System.out.println(" ---- ");
			return -1;
		}
	
		return best;
	}
	
	public double getAvgNCP(int currentLen, int tuple, int b){
		double cost = 0.0;
		if (currentLen == 0){
			return 0.0;
		}
		for(int j=0; j<currentLen; j++){
			for(int a=0; a<(dims-1); a++) {
				//cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
				if (a==0 || a==2){
					cost += (double)Math.abs(buckets[b][j][a] - map[tuple][a]) / (double)(cardinalities[a]-1);
				}else{
					if (buckets[b][j][a] != map[tuple][a])
						cost+=(double)(1.0)/(double)(cardinalities[a]-1);
				}
			}
		}
		//double denom = (double)((bucketSize - empty[b]) * canFit * (dims-1));
		double denom = (double)(bucketSize * (dims-1));
		return cost/denom;
	}
	
	public double getAvgNCP(int[] empty, ArrayList<Integer> tupSet, int offset, int b){
		double cost = 0.0;
		int canFit = tupSet.size() - offset;
		if (empty[b] < canFit) {
			canFit = empty[b];
		}
		for(int i=offset; i<(offset+canFit); i++){
			for(int j=0; j<(bucketSize - empty[b]); j++){
				for(int a=0; a<(dims-1); a++) {
					//cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					if (a==0 || a==2){
						cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a]) / (double)(cardinalities[a]-1);
					}else{
						if (buckets[b][j][a] != map[tupSet.get(i)][a])
							cost+=(double)(1.0)/(double)(cardinalities[a]-1);
					}
				}
			}
		}
		double denom = (double)((bucketSize - empty[b]) * canFit * (dims-1));
		return cost/denom;
	}

	public double getMaxNCP(int[] empty, ArrayList<Integer> tupSet, int offset, int b){
		double cost = 0.0;
		double maxcost = 0.0;
		int canFit = tupSet.size() - offset;
		if (empty[b] < canFit) {
			canFit = empty[b];
		}
		for(int i=offset; i<offset+canFit; i++){
			for(int j=0; j<(bucketSize - empty[b]); j++){
				cost = 0.0;
				for(int a=0; a<(dims-1); a++) {
					//cost = (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					if (a==0 || a==2){
						cost += (double)Math.abs(buckets[b][j][a] - map[tupSet.get(i)][a])/(double)(cardinalities[a]-1);
					}else{
						if (buckets[b][j][a] != map[tupSet.get(i)][a])
							cost += (1.0)/(double)(cardinalities[a] - 1.0);
					}
				}
				if (cost > maxcost)
					maxcost = cost;
			}
		}
		return maxcost/((double)(dims-1));
	}

	public void printBuckets(){
		// Print tuples in buckets:
		for(int b=0; b<buckNum; b++){
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
			System.err.println("Error: "+SA+" is out of range [0-"+dims+"]");
			//SA = dims;
			System.exit(-1);
		}

		//count frequencies:
		System.out.println("2. count frequencies");
		sTuples = new HashMap<Short, ArrayList<Integer>>();
		int maxFreq = 1;
		for (int i = 0; i < tuples;i++) {
			short value = map[i][SA];
			if (sTuples.containsKey(value)) {
				int tempFreq = sTuples.get(value).size() + 1;
				sTuples.get(value).add(i);
				if (maxFreq < tempFreq){
					maxFreq = tempFreq;
				}
			}else{
				ArrayList<Integer> tl = new ArrayList<Integer>();
				tl.add(i);
				sTuples.put(value, tl);
			}
		}
		
		int frSize = sTuples.size();
		//put SA values in an array to sort them!
		sArray = new short[frSize];
		int ii = 0;
		for (Map.Entry<Short, ArrayList<Integer>> entry : sTuples.entrySet()) {
			sArray[ii++] = entry.getKey();
		}
		
		FileWriter fw = null;
		try{
			fw = new FileWriter("./stats_"+inFile+".txt",false); //true == append
			fw.write("SA\t"+"freq\n");
			for (Map.Entry<Short, ArrayList<Integer>> entry : sTuples.entrySet()) {
				fw.write(entry.getKey()+"\t"+entry.getValue().size()+"\n");
			}
		}catch(Exception e){
			System.err.println("Error when printing SA stats: "+e.getMessage());
		}finally{
			try{
				if(fw != null) fw.close();
			}catch(Exception e){
				//ignore.
			}
		}
		
		//Sort SA values wrt frequency:
		quickSortSA(0, frSize-1); //sArray: from least to most frequent.
		
		//eligibility:
		System.out.println("4. Find bucket capacity (c)");
		if (b_param == 0){
			
			// Let c = GCD:
			bucketSize = gcd_all();
			System.out.println("Bucket Size: c = "+bucketSize);
			// if GCD is 1: bad utility.
			if (bucketSize==1){
				System.out.println("0-Likeness can only be achieved with n buckets of size 1 (c=1, |B|=n).");
			}
			//Create buckets:
			buckNum = (tuples / bucketSize); // div
			if (tuples % bucketSize > 0){ // mod
				buckNum++; // 1 more bucket.
			}
			buckets = new short[buckNum][bucketSize][dims];
			//Place tuples in the buckets!
			int freq = 0;
			int pos = 0; //position in the bucket
			int bi = 0; //current bucket index
			for (int i=0; i<frSize; i++){ //parse sArray
				ArrayList<Integer> tmp = sTuples.get(sArray[i]);
				freq = tmp.size();
				for (int j=0; j<freq; j++){
					buckets[bi][pos] = map[tmp.get(j)];
					pos++;
					if(pos == bucketSize){
						bi++; // next bucket
						pos = 0;
					}
				}
			}
			//printBuckets(); //For debugging.
			
			System.out.println("7. return buckets");
			return buckets;
		}
		
		//Else: b_param != 0 :

		// maximim c is maxFreq:
		int maxC = tuples; // maxFreq; //It is <=n. <<=========== NOTE: I CHANGED 'maxFreq' TO 'tuples' (WICH IS 'n', i.e., the DB size). This is only to check what happens for VERY LARGE BETA values!!!! 08/july/2019
		// minimun c is GCD:
		int minC = gcd_all(); //It is >=1.
		bucketSize = minC;
		//Find a Bucket Size c, such that: maxBeta <= b_param:
		System.out.println("minC="+minC+" maxC="+maxC);
		for(int c = maxC; c >= minC; c--){
			//System.out.println("c="+c);
			double maxRatio = 0.0; //maxRatio=max{|Bi|/ni}.
			double ratio = 0.0; //temp variable
			int part = 0; //temp variable
			int pos = 0; //temp position in bucket.
			int ni = 0;
			for (int i=0; i<frSize; i++){ //parse sTuples from sArray
				ni = sTuples.get(sArray[i]).size();
				if ((ni % c) == 0){
					ratio = 1.0/((double)c); //Because |Bi|=ni/c => (ni/c)/ni = 1/c.
				}else{
					if (ni < (c-pos)){ //fits
						pos = pos + ni;
						ratio = 1.0/((double)ni); //Because |Bi|=1.
					}else{
						part = ni - (c-pos);
						pos = (part % c);
						if (pos == 0){
							ratio = ((double)(1+(ni/c))) / ((double)ni); // (1+(ni/c))/ni = (ni+c)/ni*c
							pos=0;
						}else{
							ratio = ((double)(2+(ni/c))) / ((double)ni); // (2+(ni/c))/ni = (ni+2c)/ni*c
						}
					}
				}
				//System.out.println("ni="+ni+" ratio="+ratio);
				if(ratio > maxRatio){
					maxRatio = ratio;
				}
			}
			double maxBeta = (((double)c) * maxRatio) - 1.0;
			//System.out.println("c="+c+" maxRatio="+maxRatio+" beta>="+maxBeta);
			if (maxBeta <= b_param){
				bucketSize = c; //Found the right bucket size!
				break;
			}
		}
		System.out.println("Bucket Size: c = "+bucketSize);
		//Create buckets:
		buckNum = (tuples / bucketSize); // div
		if (tuples % bucketSize > 0){ // mod
			buckNum++; // 1 more bucket.
		}
		buckets = new short[buckNum][bucketSize][dims];
		//Place tuples in the buckets!
		int ni = 0; //support(si)
		int pos = 0; //position in the bucket
		int bi = 0; //current bucket index
		if (bucketSize == 1){
			for (int i=0; i<frSize; i++){ //parse sArray
				ArrayList<Integer> tmp = sTuples.get(sArray[i]);
				for (int j=0; j<(tmp.size()); j++){
					buckets[bi][0] = map[tmp.get(j)];
					bi++;
				}
			}
		}else{ //bucketSize > 1:
			ArrayList<Integer> tmp;
			for (int i=0; i<frSize; i++){ //First SA that fit exactly!
				tmp = sTuples.get(sArray[i]);
				ni = tmp.size();
				if ((ni % bucketSize) == 0) {
					while (!(tmp.isEmpty())) {
						if (pos == 0){
							buckets[bi][0] = map[tmp.get(0)];
							pos++;
							tmp.remove(0);
							if (tmp.isEmpty()){
								break;
								
							}
						}
						//select best record:
						int rec=0; double cost = 0;
						for (int j=0; j<tmp.size(); j++){
							double ccost = getAvgNCP(pos, tmp.get(j), bi);
							if (ccost > cost){
								rec = j;
								cost = ccost;
							}
						}
						buckets[bi][pos] = map[tmp.get(rec)];
						pos++;
						tmp.remove(rec);
						if(pos == bucketSize){
							bi++; // next bucket
							pos = 0;
						}
					}
				}
			}
			for (int i=0; i<frSize; i++){ //Then, SA that don't fit exactly!
				tmp = sTuples.get(sArray[i]);
				ni = tmp.size();
				if ((ni % bucketSize) != 0){
					while (!(tmp.isEmpty())) {
						if (pos == 0){
							buckets[bi][0] = map[tmp.get(0)];
							pos++;
							tmp.remove(0);
							if (tmp.isEmpty()){
								break;
								
							}
						}
						//select best record:
						int rec=0; double cost = 0;
						for (int j=0; j<tmp.size(); j++){
							double ccost = getAvgNCP(pos, tmp.get(j), bi);
							if (ccost > cost){
								rec = j;
								cost = ccost;
							}
						}
						buckets[bi][pos] = map[tmp.get(rec)];
						pos++;
						tmp.remove(rec);
						if(pos == bucketSize){
							bi++; // next bucket
							pos = 0;
						}
					}
				}
			}
			//Dummy Tuples:
			if((pos<bucketSize)&&(bi<buckNum)){ //Last bucket is half-full
				for(int i=pos; i<bucketSize; i++){
					for(int j=0; j<dims; j++){
						if (j == SA){
							buckets[bi][i][j] = buckets[bi][pos-1][SA]; //same SA as majority in bucket.
						}else{
							buckets[bi][i][j] = -1; //todo: ignore it while generalizing!
						}
					}
					//new dummy tuple.
				}
			}
		}
		
		//printBuckets(); //For debugging.
		
		//sTuples.clear();
		
		
		System.out.println("7. return buckets");
		return buckets;
	}
}