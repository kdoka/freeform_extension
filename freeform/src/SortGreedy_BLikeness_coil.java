import java.io.IOException;
import java.util.*;
import java.io.FileWriter;


/*includes the reading of tuples and their assignment into buckets*/
public class SortGreedy_BLikeness_coil {

	private static final double BIG = Double.MAX_VALUE;
	private static final boolean RANGE = false;
	private static final boolean MIXED = true;
	//private static final boolean OPTIMIZATION = true;
	private static double maxCost = BIG;//0.5436267458828434;
	//static byte[] cardinalities = {79, 2, 17, 6, 9, 10, 83, 51};//census
	static byte[] cardinalities ={10,6,6,10,41};//Coil2000
	//***b*like***//
	//static int k;// = 10;
	static double b_param;// = 10;
	static int SA;// 0 - 7.
	//***b*like***//
	static int dims = 5; //3
	static int tuples;// = 10000;
	static int origTuples;// = tuples;
	static int bucket_size;//c;
	static int buckNum;//tuples/c
	static int partition_size;//size of bucket partitions.
	static int partition_function;//type of bucket partitioning.
	static boolean q;
	static int parts;//number of partitions per bucket.
	//static int offset=0;
	static byte[] dimension = new byte[dims];
	static short[][] map;// = new short[tuples][dims];
	static short[][][] buckets; // = new short[tuples/c][c][dims];
	static LinkedList<Integer>[][] distinctValuesPerAssign;
	static LinkedList<Integer>[] distinctValues1;
	// = (LinkedList<Short>[][]) new LinkedList[chunk_size][dims];
	static int[][][] MinMaxPerAssign;
	static int[][] MinMaxPerAttribute; //needed for queries.
	static int[][] final_assignment;// = new int[tuples][k];
	static LinkedList<Integer> chunk_sizes = new LinkedList<Integer>();
	static HeapNode[] edges;
	static int edge_size;
	static double threshold;

	//*******************************************//
	//METHODS THAT PERFORM ARRAY-PROCESSING TASKS//
	//*******************************************//

	private static int[] greedyAssign(double[][] array, int[] assignment, int chunk_size) {

		int[] jToi = new int[chunk_size];
		Arrays.fill(assignment, -1);
		Arrays.fill(jToi, -1);
		for (int i=0; i<edge_size; i++){
			HeapNode node = (HeapNode) edges[i];
			if (assignment[node.getI()]==-1 && jToi[node.getJ()]==-1){
				assignment[node.getI()]=node.getJ();
				jToi[node.getJ()]=node.getI();
			}
		}
		//long backtrace_start = System.currentTimeMillis();
		for (int index = 0; index<chunk_size; index++){
			boolean swap_done = false;
			if (assignment[index]==-1){
				System.out.println(index);
				for (int step=1; step<chunk_size/2; step++){
					int possible_i = (chunk_size+index-step)%chunk_size;
					int possible_j = assignment[possible_i];
					if (possible_j != -1){
						if (array[index][possible_j]!=BIG){
							for (int l=0; l<chunk_size; l++){
								if (jToi[l]==-1 && array[possible_i][l]!=BIG){
									assignment[index]=possible_j;
									assignment[possible_i]=l;
									jToi[l]=possible_i;
									jToi[possible_j]=index;
									swap_done = true;
									break;
								}
							}
							/*for (int step2 = 1; step2<chunk_size/2; step2++){
							 int l = (chunk_size+possible_j-step2)%chunk_size;

							 if (jToi[l]==-1 && array[possible_i][l]!=BIG){
							 assignment[index]=possible_j;
							 assignment[possible_i]=l;
							 jToi[l]=possible_i;
							 jToi[possible_j]=index;
							 swap_done = true;
							 break;
							 }

							 l=(possible_j+step2)%chunk_size;
							 if (jToi[l]==-1 && array[possible_i][l]!=BIG){
							 assignment[index]=possible_j;
							 assignment[possible_i]=l;
							 jToi[l]=possible_i;
							 jToi[possible_j]=index;
							 swap_done = true;
							 break;
							 }
							 }	*/
						}
					}
					if (swap_done)
						break;
					else {
						possible_i = (index+step)%chunk_size;
						possible_j = assignment[possible_i];
						if (possible_j != -1){
							if (array[index][possible_j]!=BIG){
								for (int l=0; l<chunk_size; l++){
									if (jToi[l]==-1 && array[possible_i][l]!=BIG){
										assignment[index]=possible_j;
										assignment[possible_i]=l;
										jToi[l]=possible_i;
										jToi[possible_j]=index;
										swap_done = true;
										break;
									}
								}
								/*for (int step2 = 1; step2<chunk_size/2; step2++){
								 int l = (chunk_size+possible_j-step2)%chunk_size;

								 if (jToi[l]==-1 && array[possible_i][l]!=BIG){
								 assignment[index]=possible_j;
								 assignment[possible_i]=l;
								 jToi[l]=possible_i;
								 jToi[possible_j]=index;
								 swap_done = true;
								 break;
								 }

								 l=(possible_j+step2)%chunk_size;
								 if (jToi[l]==-1 && array[possible_i][l]!=BIG){
								 assignment[index]=possible_j;
								 assignment[possible_i]=l;
								 jToi[l]=possible_i;
								 jToi[possible_j]=index;
								 swap_done = true;
								 break;
								 }
								 }*/
							}
						}
					}
					if (swap_done)
						break;
				}
			}
		}

		return assignment;
	}

	//*******************************************//
	//	METHODS to pre-process the data			//
	//*******************************************//

	public static void dimension_sort() { // sort the dimensions according to their effect to GCP
		for (int i = 0;i < dims-1;i++) {
			dimension[i] = (byte) i;
		}
		boolean swapped = true;
		int i=0;
		byte temp;
		while (swapped) {
			swapped = false;
			i++;
			for (int j = 0; j < dimension.length-i;j++) { // smaller range, put it the front
				if (cardinalities[dimension[j]] > cardinalities[dimension[j+1]]) {
					temp = dimension[j];
					dimension[j] = dimension[j+1];
					dimension[j+1] = temp;
					swapped = true;
				}
			}
		}
	}

	// return true if data x > data y
	public static boolean compare_data(short[] x, short[] y) {
		for (int i = 0;i < dims;i++) {
			if (i != dims-1){ //compare w.r.t. all except SA!!!
				if (x[dimension[i]] > y[dimension[i]]) {
					/*if (map[x][dimension[0]]<map[y][dimension[0]])
					 System.out.println("oops");*/
					return true;
				} else if (x[dimension[i]] < y[dimension[i]]) {
					/*if (map[x][dimension[0]]>map[y][dimension[0]])
					 System.out.println("oops");*/
					return false;
				}
			}
		}
		/*if (map[x][dimension[0]]>map[y][dimension[0]])
			System.out.println("oops");*/
		return false;
	}

	// sort edges[] by lexicographical order
	public static void qSort(int left, int right){
		int i = left, j = right;
		int ref = left + (right-left)/2;
		HeapNode pivot = edges[ref];
		HeapNode temp2 ;
		while (i <= j) {
			while (compare(pivot, edges[i]))
				i++;
			while (compare(edges[j], pivot))
				j--;
			if (i <= j) {
				temp2=edges[i];
				edges[i]=edges[j];
				edges[j]=temp2;

				i++;
				j--;
			}
		};
		// recursion
		if (left < j)
			qSort(left, j);
		if (i < right) {
			qSort(i, right);
		}

	}

	public static boolean compare(Object o1, Object o2) {
		if (((HeapNode) o1).getCost()>((HeapNode)o2).getCost())
			return true;
		else
			return false;
	}

	public static void quickSortArray(int left, int right, double[] a) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		double pivot = a[ref];
		double temp2;
		while (i <= j) {
			while (pivot > a[i])
				i++;
			while (a[j] > pivot)
				j--;
			if (i <= j) {
				temp2=a[i];
				a[i]=a[j];
				a[j]=temp2;

				i++;
				j--;
			}
		};
		// recursion
		if (left < j)
			quickSortArray(left, j, a);
		if (i < right) {
			quickSortArray(i, right, a);
		}
	}
	// sort data by lexicographical order
	public static void quickSortBucket(int left, int right, int b) {
		int i = left, j = right;
		int ref = left + (right-left)/2; //i.e., (right+left)/2
		short[] pivot = buckets[b][ref];
		short[] temp2 = new short[dims];
		while (i <= j) {        
			while (compare_data(pivot, buckets[b][i]))
				i++;
			while (compare_data(buckets[b][j], pivot))
				j--;
			if (i <= j) {
				temp2=buckets[b][i];
				buckets[b][i]=buckets[b][j];
				buckets[b][j]=temp2;

				i++;
				j--;
			}
		};    
		// recursion
		if (left < j)
			quickSortBucket(left, j, b);
		if (i < right) {
			quickSortBucket(i, right, b);
		}
	}

	/*
	 * Option 0: bucket_partition.
	 * Partitions each bucket, while maintaining SA distributions
	 * Note: This is the correct one!
	 */
	static void bucket_partition(int b_size, LikenessBuckets bk){
		int parts = b_size/partition_size; //#partitions per bucket
		double ratio; int offsetB; int offsetP;
		int chunk_size; int partSA;
		short[][] tmpBucket;
		int[] loaded;

		if(parts == 0){
			parts = 1;
			chunk_sizes.add(b_size);
		}else{
			for(int i=0; i<parts-1; i++){
				chunk_sizes.add(partition_size);
			}
			if ((b_size % partition_size) <= (partition_size/2)){
				chunk_sizes.addFirst(partition_size + (b_size % partition_size));
			}else{
				parts++;
				chunk_sizes.add(partition_size);
				chunk_sizes.addLast(b_size % partition_size);
			}
		}
		loaded = new int[parts];

		for(int b=0; b<buckNum; b++){//for every bucket
			tmpBucket = new short[bucket_size][dims];
			ArrayList<Integer> chg = bk.changeSA.get(b);
			int first = 0; int last;
			for(int i=0; i<parts; i++){
				loaded[i] = 0;
			}
			for(int i=0; i<chg.size(); i++){//for every SA
				last = chg.get(i);
				if((last-first) > 1){
					//System.out.println("sorting bucket["+b+"]["
					//				   +first+"--"+last+"]");
					quickSortBucket(first, last, b);
				}
				offsetB = first;
				offsetP = 0;
				for(int jj=0; jj<parts; jj++){//for every partition
					chunk_size = chunk_sizes.get(jj);
					//ratio = ((double)chunk_size) / ((double)bucket_size);
					//partSA = (int)(Math.ceil(ratio * ((float)(last-offsetB+1))));
					partSA = (int)Math.ceil(((double)(last-offsetB+1)) / ((double)(parts-jj)));
					if (partSA > (chunk_size - loaded[jj]))
						partSA = chunk_size - loaded[jj];
					for(int j=0; j<partSA; j++){
						tmpBucket[offsetP + loaded[jj] + j] = buckets[b][offsetB + j];
					}
					loaded[jj] += partSA;
					offsetB += partSA;
					offsetP += chunk_size;
				}
				first = last+1;
			}
			buckets[b] = tmpBucket;
		}
	}

	/*
	 * Option 1: bucket_partition2.
	 * Partitions each bucket, maintaining the distribution of the most freq SA
	 * while partitioning the rest  based only on their QIDs, ignoring SAs.
	 * Note: It preserves better data utility, but may lead to a deadlock
	 * (chunk assignments that have 50% the same SAs.)
	 */
	static void bucket_partition2(int b_size, LikenessBuckets bk){
		int parts = b_size/partition_size; //#partitions per bucket
		float ratio; int offsetB; int offsetB2; int offsetP;
		int chunk_size; int partSA; int partSA2;
		short[][] tmpBucket;

		if(parts == 0){
			parts = 1;
			chunk_sizes.add(b_size);
		}else{
			for(int i=0; i<parts-1; i++){
				chunk_sizes.add(partition_size);
			}
			if ((b_size % partition_size) <= (partition_size/2)){
				chunk_sizes.addFirst(partition_size + (b_size % partition_size));
			}else{
				parts++;
				chunk_sizes.add(partition_size);
				chunk_sizes.addLast(b_size % partition_size);
			}
		}

		for(int b=0; b<buckNum; b++){//for every bucket
			tmpBucket = new short[bucket_size][dims];
			//tuples of most freq SA:
			ArrayList<Integer> chg = bk.changeSA.get(b);
			int first = 0; int last = chg.get(0);
			if((last-first) > 0){ // (last-first+1) > 1.
				//System.out.println("sorting bucket["+b+"]["
				//				   +first+"--"+last+"]");
				quickSortBucket(first, last, b);
			}
			//remaining tuples:
			if((bucket_size - last) > 3){ // (bucket_size-1 -last-1) > 1.
				//System.out.println("sorting bucket["+b+"]["
				//				   +(last+1)+"--"+(bucket_size-1)+"]");
				quickSortBucket(last+1, bucket_size-1, b);
			}
			offsetB = 0;
			offsetP = 0;
			offsetB2 = last+1;
			for(int jj=0; jj<parts; jj++){//for every partition
				chunk_size = chunk_sizes.get(jj);
				//ratio = ((float)chunk_size) / ((float)bucket_size);
				//partSA = Math.round(ratio * ((float)(last-first+1))); //freq SA
				partSA = (int)Math.ceil(((double)(last-offsetB+1)) / ((double)(parts-jj)));
				partSA2 = chunk_size - partSA; //remaining SAs
				for(int j=0; j<partSA; j++){
					tmpBucket[offsetP + j] = buckets[b][offsetB + j];
				}
				offsetB += partSA;
				for(int j=0; j<partSA2; j++){
					tmpBucket[offsetP + partSA + j] = buckets[b][offsetB2 + j];
				}
				offsetB2 += partSA2;
				offsetP += chunk_size;
			}
			buckets[b] = tmpBucket;
		}
	}

	static double rangeQueries(double s, int lambda, double[][] rs, int times, double[] errArray){
		double error = 0.0;
		double absError=0.0;
		double genMx, genMn;
		double sele, range, min, max, random, randomPos, r, overlap;
		int[] distVals = new int[buckNum];
		ArrayList<Integer> attrOptions = new ArrayList<Integer>();
		ArrayList<Integer> choices = new ArrayList<Integer>();
		//static final int[] attrNumber = new int[]{0, 1, 2, 3, 4, 5 , 6};

		sele = Math.pow(s, 1.0/((double)lambda+1.0));

		System.out.println("lambda ="+lambda+": ");

		for (int tm=0; tm<times; tm++){
			for (int qi=0; qi<dims-1; qi++){ //INITIALIZATION:
				rs[qi][0] = (double)MinMaxPerAttribute[qi][0]-1;//<min SA value.
				rs[qi][1] = (double)MinMaxPerAttribute[qi][1]+1;//>max SA value
			}
			//SA:
			range = (double)(MinMaxPerAttribute[SA][1] - MinMaxPerAttribute[SA][0]) * sele;
			min = (double)MinMaxPerAttribute[SA][0];//min SA value.
			max = (double)MinMaxPerAttribute[SA][1] - range;//max SA value - range.
			random = new Random().nextDouble();
			rs[dims-1][0] = min + (random * (max-min));
			rs[dims-1][1] = rs[dims-1][0] + range;
			//others:
			attrOptions.clear();
			choices.clear();
			for (int i=0; i<dims-1; i++){ attrOptions.add(i); } //originally

			int temporary = dims-1;
			for (int ch=0; ch<lambda; ch++){
				int next = new Random().nextInt(temporary);
				//System.out.println(next+" : "+(attrOptions.get(next)));
				choices.add(attrOptions.get(next));
				attrOptions.remove(next);
				temporary--;
			}
			/*
			System.out.print("Keeping:= ");
			for (int i=0; i<choices.size(); i++){
				System.out.print(choices.get(i)+" ");
			}
			System.out.print("-- Discarding: ");
			for (int i=0; i<attrOptions.size(); i++){
				System.out.print(attrOptions.get(i)+" ");
			}
			if (choices.size() != lambda){
				System.out.println("WTF: "+choices.size()+" "+lambda);
				System.exit(0);
			}
			System.out.println();
			 */
			for (int ii=0; ii<lambda; ii++){
				int i = (int)choices.get(ii);
				//Categorical:
				range = (double)Math.round((double)(cardinalities[i]) * sele);
				if (range < 1.0) range = 1; //select at least 1 value.
				min = (double)MinMaxPerAttribute[i][0];//min attribute value.
				//max = (double)MinMaxPerAttribute[i][1] - range;//max attribute value - range.
				random = new Random().nextDouble();
				randomPos = (double)Math.round(random * ((double)(cardinalities[i]) * (1.0-sele)));
				rs[i][0] = min + randomPos;
				rs[i][1] = rs[i][0] + range;

			}
			/*for (int i=0; i<lambda; i++)
			 System.out.println(i+"s="+s+"sele="+sele+":["+MinMaxPerAttribute[i][0]+","
			 +MinMaxPerAttribute[i][1]+"]"+rs[i][0]+" "+rs[i][1]);
			 System.out.println(" ");*/
			double cnt = 0; double anonCnt = 0;
			for (int j=0; j<final_assignment.length; j++){ //not origTuples!
				r = 1.0; overlap=1.0;
				boolean inRange = true;
				boolean genRange = true;
				for (int index=0; index<lambda; index++){
					int i = (int)choices.get(index);
					genMx = indexToTupleMapping(final_assignment[j][0])[i];
					genMn = indexToTupleMapping(final_assignment[j][0])[i];

					for (int bi=0; bi<buckNum; bi++){
						//Generalization Range:
						if (final_assignment[j][bi] < origTuples){ //not a dummy:
							distVals[bi] = indexToTupleMapping(final_assignment[j][bi])[i];
							if (genMx < distVals[bi])
								genMx = distVals[bi];
							if ((genMn > distVals[bi]) || (genMn==-1))
								if (-1 != distVals[bi])
									genMn = distVals[bi];
						}
					}
					if ((genMx<rs[i][0])||(genMn>rs[i][1])){//gen out of query range.
						//if (!((genMx<rs[i][1])&&(genMn>rs[i][0]))){//gen not within range.
						genRange = false;
						//break;
					}else{ //there is some overlap:

						if ((rs[i][0]<=genMn)&&(genMx<=rs[i][1])){
							overlap = 1.0;
						}else{
							int tempCnt = 0;
							for (int it=0; it<distVals.length; it++){
								if((distVals[it]>=rs[i][0])&&(distVals[it]<=rs[i][1]))
									tempCnt++;
							}
							overlap = ((double)tempCnt) / ((double)distVals.length);
						}

						r = r * overlap;
						//System.out.println("r="+r);
					}

					if ((indexToTupleMapping(final_assignment[j][0])[i]<rs[i][0])||
							(indexToTupleMapping(final_assignment[j][0])[i]>rs[i][1])){
						//orig out of range:
						inRange = false;
						break;
					}
				}

				if ((indexToTupleMapping(final_assignment[j][0])[SA]<rs[SA][0])||
						(indexToTupleMapping(final_assignment[j][0])[SA]>rs[SA][1])){//out of range
					inRange = false;
					genRange = false;
				}

				if (true == inRange){
					//System.out.println("orig +1 ");
					cnt++;
				}
				if (true == genRange){
					//System.out.println("anon +r "+r);
					anonCnt += r;
				}
			}
			//System.out.println(cnt+" "+anonCnt);
			if (cnt!=0){
				//System.out.print(tm+" ");
				error += (((double)Math.abs(anonCnt - cnt)) / ((double)cnt));
				absError += (((double)Math.abs(anonCnt - cnt)));
				errArray[tm] = (((double)Math.abs(anonCnt - cnt)) / ((double)cnt));
			}else{
				tm--; //repeat.
			}

		}
		quickSortArray(0, times-1, errArray);
		double median = (errArray[(times/2)-1]);
		System.out.println("query error (sel="+s+", lambda="+lambda+"): mean rel error="
				+((error/(double)times))+" abs error="+(absError/times)
				+" median rel error="+median);
		System.out.println("min="+errArray[0]+" max="+errArray[times-1]);
		return median;

		/*
		 double mean = 0.0;
		 for (i=0; i<errArray.size(); i++){
		 mean += errArray[i];
		 }
		 mean = mean / errArray.size();
		 return mean;
		 */

	}


	//***********//
	//MAIN METHOD//
	//***********//

	public static void main(String[] args) 	{

		if (args.length!=7){
			System.out.println("\nUsage:   java SortGreedy_BLikeness inFile n SA beta part_size part_option q");
			System.out.println("\t inFile: input file name (path included).");
			System.out.println("\t n: number of tuples in inFile.");
			//System.out.println("\t d: dimensionality of the dataset.");
			System.out.println("\t SA: index of sensitive attribute [0 -- d-1].");
			System.out.println("\t l: beta: B-likeness parameter.");
			System.out.println("\t part_size: size of the bucket partitions.");
			System.out.println("\t part_option: 0 (safer, keeps all SAs distributions), or");
			System.out.println("\t              1 (better utility, but may cause problems), or ");
			System.out.println("\t              2 (no bucket partitioning).\n");
			System.out.println("\t q: 1 with queries, 0 without\n");
			//			System.out.println("\t th: distance threshold to place chunk in bucket, in [0, 1].");
			return;
		}

		String inputFile = args[0];
		tuples = Integer.parseInt(args[1]);  // n
		//dims = Integer.parseInt(args[2]); //d
		SA = Integer.parseInt(args[2]); //Sensitive Attribute (0 - 7).
		b_param = Double.parseDouble(args[3]); // beta
		partition_size = Integer.parseInt(args[4]);
		partition_function = Integer.parseInt(args[5]);
		q = Boolean.parseBoolean(args[6]);
		//		threshold = Double.parseDouble(args[7]);
		origTuples = tuples;

		/*
		int modl = (tuples % l_param);
		if (modl > 0){
			//change n (#tuples), so that it is divided by l:
			tuples = tuples + (l_param - modl);
			for (int i=0; i<dims-1; i++)
				cardinalities[1]++;
		}
		map = new short[tuples][dims];
		buckets = new short[l_param][tuples/l_param][dims];
		 */
		map = new short[tuples][dims];
		MinMaxPerAttribute = new int[dims][2];

		long startTime = System.currentTimeMillis();
		try {
			CensusParser tp = new CensusParser(inputFile, dims);
			int i=0;
			while (tp.hasNext()){
				map[i++]=tp.nextTuple2();
				for (int j=0; j<dims; j++){
					if (map[i-1][j] < MinMaxPerAttribute[j][0]){ //min
						MinMaxPerAttribute[j][0] = map[i-1][j];
					}
					if (map[i-1][j] > MinMaxPerAttribute[j][1]){ //max
						MinMaxPerAttribute[j][1] = map[i-1][j];
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		/*
		if (modl > 0){
			//add dummy tuples:
			for(int i=(tuples-(l_param-modl)); i<tuples; i++){
				for(int j=0; j<dims; j++){
					if (j == dims-1){
						map[i][j] = -1; //unique SA
					}else{
						map[i][j] = 0;
					}
				}
			}
		}
		 */
		long midTime = System.currentTimeMillis();
		dimension_sort();//sort the dimensions
		LikenessBuckets bk = new LikenessBuckets(b_param, tuples, dims, map, buckets, 0, inputFile);
		buckets = bk.bucketization(SA);
		//bk.printBuckets();

		long bucketEndTime = System.currentTimeMillis();
		System.out.println("Time of reading dataset: "+(midTime - startTime)+" miliseconds.");
		System.out.println("Time of creating buckets: "+(bucketEndTime - midTime)+" miliseconds.");
		//bk.printBuckets();

		map=null; //delete map
		System.gc();


		bucket_size = bk.bucketSize;//bucket capacity (c).
		buckNum = bk.buckNum; //number of buckets (|B|).
		System.out.println("Number of buckets:"+buckNum);
		//update "tuples" number, taking into account the dummies:
		tuples = (bucket_size * buckNum);

		final_assignment = new int[tuples][buckNum];

		double distortion = 0.0;
		int chunk_size;

		/*
		 * Sort groups of same-SA tuples in each bucket, wrt QIDs.
		 * Then form bucket partitions, keeping SA distributions.
		 */
		if (partition_function == 0)
			bucket_partition(bucket_size, bk); //keep all SAs distrubutions.
		else if (partition_function == 1){
			bucket_partition2(bucket_size, bk); //only keep 1st SA distribution.
		} //else NO_PARTITION //default.
		System.gc();
		//bk.printBuckets();

		if((partition_function == 0) || (partition_function == 1)){ //partitioned buckets:

			for (int bucket_index=0; bucket_index<buckNum; bucket_index++){

				int chunk_offset = 0;
				for (int chunk_index=0; chunk_index<chunk_sizes.size(); chunk_index++){
					chunk_size = chunk_sizes.get(chunk_index);

					//edges = new HeapNode[chunk_size*chunk_size*2];
					edges = new HeapNode[chunk_size*chunk_size];

					//we need SA, too!
					if (MIXED){
						MinMaxPerAssign = new int[chunk_size][dims-1][2];
						distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];
					}else{
						if (RANGE){
							MinMaxPerAssign = new int[chunk_size][dims-1][2];
							distinctValues1 = (LinkedList<Integer>[]) new LinkedList[chunk_size];
						}else{
							distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];
						}
					}
					//HeapComparator hc = new HeapComparator();
					double[][] array = computeCostMatrix(buckets[bucket_index],buckets[(bucket_index+1)%buckNum], bucket_index*bucket_size, chunk_offset, chunk_size);
					int[] assignment = new int[array.length];
					int times = 0;

					while (++times<buckNum){
						//qSort(0, chunk_size*chunk_size-1);
						qSort(0, edge_size-1);
						greedyAssign(array, assignment, chunk_size);//Call SortGreedy algorithm.

						for (int i=0; i<assignment.length; i++){
							final_assignment[i+chunk_offset+bucket_index*bucket_size][times] = bucketToIndexMapping((bucket_index+times)%buckNum,(chunk_offset+assignment[i]));
							if (MIXED){
								findSet_mixed(i, buckets[(bucket_index+times)%buckNum][chunk_offset+assignment[i]] );
							}else{
								if (RANGE)
									findSet_numerical(i, buckets[(bucket_index+times)%buckNum][chunk_offset+assignment[i]]);
								else
									findSet(i,buckets[(bucket_index+times)%buckNum][chunk_offset+assignment[i]]);
							}
						}
						if (times!=buckNum-1)
							recomputeCostMatrix(array, (bucket_index+times+1)%buckNum, chunk_offset, chunk_size);
					}
					for (int i=0; i<chunk_size; i++){
						if (MIXED){
							distortion += NCP_mixed(MinMaxPerAssign[i], distinctValuesPerAssign[i]);
						}else
							if (RANGE)
								distortion += NCP_numerical(MinMaxPerAssign[i]);
							else
								distortion += NCP(distinctValuesPerAssign[i]);
					}
					chunk_offset += chunk_size;
				}
			}
		}else{ //No partitioning:
			for (int bucket_index=0; bucket_index<buckNum; bucket_index++){

				//edges = new HeapNode[bucket_size*bucket_size*2];
				edges = new HeapNode[bucket_size*bucket_size];
				//we need SA, too!
				if (MIXED){
					MinMaxPerAssign = new int[bucket_size][dims-1][2];
					distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[bucket_size][dims];
				}else{
					if (RANGE){
						MinMaxPerAssign = new int[bucket_size][dims-1][2];
						distinctValues1 = (LinkedList<Integer>[]) new LinkedList[bucket_size];
					}else{
						distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[bucket_size][dims];
					}
				}
				//HeapComparator hc = new HeapComparator();
				double[][] array = computeCostMatrix(buckets[bucket_index],buckets[(bucket_index+1)%buckNum], bucket_index*bucket_size, 0, bucket_size);

				int[] assignment = new int[array.length];
				int times = 0;

				while (++times<buckNum){
					//qSort(0, bucket_size*bucket_size-1);
					qSort(0, edge_size-1);
					greedyAssign(array, assignment, bucket_size);//Call SortGreedy algorithm.

					for (int i=0; i<assignment.length; i++){
						final_assignment[i+bucket_index*bucket_size][times] = bucketToIndexMapping((bucket_index+times)%buckNum, assignment[i]);
						if (MIXED){
							findSet_mixed(i, buckets[(bucket_index+times)%buckNum][assignment[i]] );
						}else{
							if (RANGE)
								findSet_numerical(i, buckets[(bucket_index+times)%buckNum][assignment[i]]);
							else
								findSet(i,buckets[(bucket_index+times)%buckNum][assignment[i]]);
						}
					}
					if (times!=buckNum-1)
						recomputeCostMatrix(array, (bucket_index+times+1)%buckNum, 0, bucket_size);
				}
				for (int i=0; i<bucket_size; i++){
					if (MIXED){
						distortion += NCP_mixed(MinMaxPerAssign[i], distinctValuesPerAssign[i]);
					}else
						if (RANGE)
							distortion += NCP_numerical(MinMaxPerAssign[i]);
						else
							distortion += NCP(distinctValuesPerAssign[i]);
				}
			}
		}//endif (partition or no_partition)

		//**** BEGIN XUE MINGQIANG **** //
		// this call returns a random assignment generated from the k-regular matching graph
		int [] rand_A = Randomization.run(final_assignment, 0, final_assignment.length, buckNum);
		//**** END XUE MINGQIANG **** //
		long endTime = System.currentTimeMillis();

		//System.out.println("The winning assignment after "+index+" runs (" + sumType + " sum) is:\n");	

		/*
		for (int i=0; i<final_assignment.length; i++){
			for (int j=0; j<buckNum; j++){
				System.out.print((final_assignment[i][j] +1)+" ");
			}
			System.out.println();
		}
		 */

		System.out.println("Time: "+(endTime - startTime)+"ms  "+"\n Distortion "+ (double)(distortion/((dims-1)*tuples)));

		System.out.println("Saving results.");
		//Save Results:
		FileWriter fw = null;
		try{
			fw = new FileWriter("./SortGreedyResults.txt",true); //true == append
			fw.write(origTuples+" "+b_param+" ");
			if((partition_function == 0) || (partition_function == 1)){
				fw.write(partition_size+" ");
			}else{
				fw.write(bucket_size+" ");
			}
			fw.write((endTime - startTime)+" "
					+((double)(distortion/((dims-1)*tuples)))+"\n");
		}catch(IOException ioe){
			System.err.println("IOException: " + ioe.getMessage());
		}finally{
			try{
				if(fw != null) fw.close();
			}catch(Exception e){
				//ignore.
			}
		}
		if (!q){
			System.out.println("Range Queries.");
			double[] selectivities = {0.05, 0.1, 0.15, 0.2, 0.25};
			double qErr = 0;
			FileWriter qw = null;
			try{
				int qtimes = 1000; //numer of random queries.
				double[] errArray = new double[qtimes];
				qw = new FileWriter("./SortGreedy_QueryError.txt",true); //true == append
				//qw.write("#tuples beta size lamda sel error\n");
				/*
			for (int i=0; i<selectivities.length; i++){
				for (int l=1; l<dims; l++){
					qw.write(origTuples+" "+b_param+" "+bucket_size+" "+
							 l+" "+selectivities[i]+" ");
					double[][] tmpres = new double[l+1][2];
					qErr = rangeQueries(selectivities[i], l, tmpres, qtimes);
					qw.write(qErr+" \n");
				}
			}
				 */
				//sel=0.1
				double[][] tmpres = new double[dims][2];

				System.out.println("Vary lambda (sel=0.1): ");
				for (int l=1; l<dims; l++){
					qw.write(origTuples+" "+b_param+" "+bucket_size+" "+
							l+" "+selectivities[1]+" ");
					for (int qi=0; qi<dims; qi++){ //INITIALIZATION:
						tmpres[qi][0] = (double)MinMaxPerAttribute[qi][0]-1;//<min SA value.
						tmpres[qi][1] = (double)MinMaxPerAttribute[qi][1]+1;//>max SA value
					}

					qErr = rangeQueries(selectivities[1], l, tmpres, qtimes, errArray);
					qw.write(qErr+" \n");
				}

				qw.write("\n");
				System.out.println("Vary selectivity (lambda=3): ");
				int l=3; //lambda = 3 first QIs.
				for (int i=0; i<selectivities.length; i++){
					qw.write(origTuples+" "+b_param+" "+bucket_size+" "+
							l+" "+selectivities[i]+" ");
					for (int qi=0; qi<dims; qi++){ //INITIALIZATION:
						tmpres[qi][0] = (double)MinMaxPerAttribute[qi][0]-1;//<min SA value.
						tmpres[qi][1] = (double)MinMaxPerAttribute[qi][1]+1;//>max SA value
					}

					qErr = rangeQueries(selectivities[i], l, tmpres, qtimes, errArray);
					qw.write(qErr+" \n");
				}

				qw.write("\n");
			}catch(IOException ioe){
				System.err.println("IOException: " + ioe.getMessage());
			}finally{
				try{
					if(qw != null) qw.close();
				}catch(Exception e){
					System.err.println(e.getMessage());
				}
			}
		}

	}
	private static double[][] computeCostMatrix(short[][] in1, short[][] in2, int offset, int first, int size) {
		edge_size = 0;

		double[][] cost = new double[size][size];
		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				double c;
				if (MIXED){
					c=NCP_mixed(in1[first+i], in2[first+j]);
				}else if (RANGE){
					c=NCP_numerical(in1[first+i], in2[first+j]);
				}else{
					c=NCP(in1[first+i], in2[first+j]);
				}
				if (c<=maxCost){
					cost[i][j]=c;
				}else{
					cost[i][j]=BIG;
				}
				HeapNode hn = new HeapNode(i,j,c);
				edges[edge_size++]=hn;

			}
			final_assignment[offset+first+i][0]= offset+first+i;
			if (MIXED){
				for (int l=0; l<dims; l++){
					LinkedList<Integer> list = new LinkedList<Integer>();
					if (in1[first+i][l] != -1){ //dummy
						list.add((int) in1[first+i][l]);
					}
					distinctValuesPerAssign[i][l] = list;
				}
			}else{
				if (RANGE){
					for (int l=0; l<dims-1; l++){
						MinMaxPerAssign[i][l][0]=in1[first+i][l];
						MinMaxPerAssign[i][l][1]=in1[first+i][l];
					}
					LinkedList<Integer> list = new LinkedList<Integer>();
					//if (in1[first+i][SA-1] != -1){ //dummy
					list.add((int) in1[first+i][dims-1]);
					//}
					distinctValues1[i] = list;
				}else{
					for (int l=0; l<dims; l++){
						LinkedList<Integer> list = new LinkedList<Integer>();
						if (in1[first+i][l] != -1){ //dummy
							list.add((int) in1[first+i][l]);
						}
						distinctValuesPerAssign[i][l] = list;
					}
				}
			}
		}
		return cost;
	}

	private static void recomputeCostMatrix(double[][] array, int bucket_index, int first, int size) {
		edge_size=0;

		for (int i=0; i<size; i++){
			for (int j=0; j<size; j++){
				if (MIXED){
					array[i][j]=NCP_mixed(buckets[bucket_index][first+j], MinMaxPerAssign[i],
							distinctValuesPerAssign[i]);
					HeapNode hn = new HeapNode(i,j,(array[i][j]));
					edges[edge_size++]=hn;
				}else if (RANGE){
					array[i][j]=NCP_numerical(buckets[bucket_index][first+j], MinMaxPerAssign[i],
							distinctValues1[i]);
					HeapNode hn = new HeapNode(i,j,(array[i][j]));
					edges[edge_size++]=hn;
				}else{
					array[i][j]=NCP(buckets[bucket_index][first+j], (distinctValuesPerAssign[i]));
					HeapNode hn = new HeapNode(i,j,(array[i][j]));
					edges[edge_size++]=hn;
				}
			}
		}
	}

	///////////////////Mixed Representation//////////////	
	private static double NCP_mixed(short[] tuple1, short[] tuple2){
		double score=0.0;

		//if (tuple1[SA] == tuple2[SA]) //Beta-likeness does not require this!.
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if((tuple1[i]==-1)||(tuple2[i]==-1)){ //Beta-likeness dummy tuple.
				return 0;
			}

			if (tuple1[i]==tuple2[i])
				score+=0;
			else 
				score+=(double)(1)/(double)(cardinalities[i]-1);

		}
		return score;
	}

	private static double NCP_mixed(short[] tuple, int[][] MinMaxPerDim, LinkedList<Integer>[] distinctValuesPerDim ){
		double score=0.0;
		int min;
		int max;

		LinkedList<Integer> distinctValues2 = distinctValuesPerDim[dims-1];
		//if (distinctValues2.contains((int)tuple[SA])) //Beta-likeness does not require this!.
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if(tuple[i]==-1){ //Beta-likeness dummy tuple.
				return 0;
			}
			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			if (!distinctValues.contains((int)tuple[i])){
				if(distinctValues.contains(-1))
					score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1); //dummy
				else
					score+=(double)(distinctValues.size())/(double)(cardinalities[i]-1);
			}else{
				if(distinctValues.contains(-1))
					score+=(double)(distinctValues.size()-2)/(double)(cardinalities[i]-1); //dummy
				else
					score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
			}

		}
		return score;
	}
	private static double NCP_mixed(int[][] MinMaxPerDim, LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;

		for (int i=0; i<dims-1; i++){

			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);

		}
		return score;
	}

	private static void findSet_mixed(int assign_number, short[] newTuple){
		LinkedList<Integer>[] dValues = distinctValuesPerAssign[assign_number];
		int[][] MinMaxPerDim = MinMaxPerAssign[assign_number];
		if(newTuple[1]==-1){ //Beta-Likeness -- not sure about this!
			return; //dummy tuple
		}
		for (int i=0; i<dims; i++){ //we need SA, too!

			LinkedList<Integer> distValuesPerDim = dValues[i];
			if (distValuesPerDim != null)
				if (!distValuesPerDim.contains((int)newTuple[i]))
					distValuesPerDim.add((int)newTuple[i]);


		}
		return;
	}


	///////////////////Range Representation//////////////			
	private static double NCP_numerical(short[] tuple1, short[] tuple2){
		double score=0.0;

		//if (tuple1[SA] == tuple2[SA]) //Beta-likeness does not require this!
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if((tuple1[i]==-1)||(tuple2[i]==-1)){ //Beta-likeness dummy tuple.
				return 0;
			}
			score+=(double)Math.abs(tuple1[i]-tuple2[i])/(double)(cardinalities[i]-1);

		}
		return score;
	}

	private static double NCP_numerical(short[] tuple, int[][] MinMaxPerDim,
			LinkedList<Integer> distinctValues2){
		double score=0.0;
		int min;
		int max;

		//if (distinctValues2.contains((int)tuple[SA]))  //Beta-likeness does not require this!
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if(tuple[i]==-1){ //Beta-likeness dummy tuple.
				return 0;
			}
			int[] distinctValues = MinMaxPerDim[i];
			min = distinctValues[0];
			max = distinctValues[1];
			if (tuple[i] < min)
				score+=(double)(max-tuple[i])/(double)(cardinalities[i]-1);
			else if ((tuple[i]>max)&&(max!=-1))
				score+=(double)(tuple[i]-min)/(double)(cardinalities[i]-1);
			else
				score+=(double)(max-min)/(double)(cardinalities[i]-1);
		}
		return score;
	}

	private static double NCP_numerical(int[][] MinMaxPerDim){
		double score=0.0;
		for (int i=0; i<dims-1; i++){
			int[] distinctValues = MinMaxPerDim[i];
			score+=(double)(distinctValues[1]-distinctValues[0])/(double)(cardinalities[i]-1);

		}
		return score;
	}

	private static void findSet_numerical(int assign_number, short[] newTuple){
		int[][] MinMaxPerDim = MinMaxPerAssign[assign_number];
		if(newTuple[0]==-1){ //Beta-Likeness -- not sure about this!
			return; //dummy tuple
		}
		for (int i=0; i<MinMaxPerDim.length; i++){
			int[] distinctValues = MinMaxPerDim[i];
			if ((newTuple[i] < distinctValues[0])||(distinctValues[0]==-1)){
				distinctValues[0]=(int) newTuple[i];
			}
			if ((newTuple[i] > distinctValues[1])||(distinctValues[1]==-1)){
				distinctValues[1]=(int) newTuple[i];
			}
		}

		return;
	}

	////////////////Set Representation///////////////
	private static double NCP(short[] tuple1, short[] tuple2){
		double score=0.0;

		//if (tuple1[SA] == tuple2[SA]) //Beta-likeness does not require this!
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if((tuple1[i]==-1)||(tuple2[i]==-1)){ //Beta-likeness dummy tuple.
				return 0;
			}
			if (tuple1[i]==tuple2[i])
				score+=0;
			else 
				score+=(double)(1)/(double)(cardinalities[i]-1);
		}
		return score;
	}

	private static double NCP(short[] tuple, LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;

		LinkedList<Integer> distinctValues2 = distinctValuesPerDim[SA];
		//if (distinctValues2.contains((int)tuple[SA])) //Beta-likeness does not require this!
		//	return BIG; //inf

		for (int i=0; i<dims-1; i++){
			if(tuple[i]==-1){ //Beta-likeness dummy tuple.
				return 0;
			}
			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			if (!distinctValues.contains((int)tuple[i])){
				if (distinctValues.contains(-1))
					score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1); //dummy
				else
					score+=(double)(distinctValues.size())/(double)(cardinalities[i]-1);
			}else{
				if (distinctValues.contains(-1))
					score+=(double)(distinctValues.size()-2)/(double)(cardinalities[i]-1); //dummy
				else
					score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
			}
		}
		return score;
	}
	private static double NCP(LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;
		for (int i=0; i<dims-1; i++){
			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			if(distinctValues.contains(-1))
				score+=(double)(distinctValues.size()-2)/(double)(cardinalities[i]-1); //dummy
			else
				score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
		}
		return score;
	}

	private static void findSet(int assign_number, short[] newTuple){
		LinkedList<Integer>[] distinctValues = distinctValuesPerAssign[assign_number];
		if(newTuple[0]==-1){ //Beta-Likeness -- not sure about this!
			return; //dummy tuple
		}
		for (int i=0; i<distinctValues.length; i++){
			LinkedList<Integer> distValuesPerDim = distinctValues[i];
			if (distValuesPerDim != null)
				if (!distValuesPerDim.contains((int)newTuple[i]))
					distValuesPerDim.add((int)newTuple[i]);
		}
		return;
	}

	private static int bucketToIndexMapping(int bucket, int i){
		return bucket*bucket_size+i;
	}

	private static short[] indexToTupleMapping(int index){
		return buckets[index/bucket_size][index%bucket_size];
	}
}
