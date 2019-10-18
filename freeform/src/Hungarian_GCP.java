import java.io.IOException;
import java.util.*;


/*if used with this hungarian change assignment in main*/
public class Hungarian_GCP {

	private static final double BIG = Double.MAX_VALUE;
	private static final boolean RANGE = false;
	private static final boolean MIXED = true;
	private static final boolean OPTIMIZATION = true;
	private static double maxCost = BIG;//0.5436267458828434;
	static byte[] cardinalities = {79, 2, 17, 6, 9, 10, 83, 51};
	//static byte[] cardinalities ={2};
	static int partition_size;// = 1000;
	static int k;// = 10;
	static int dims = 8;
	static int tuples;// = 10000;
	static int chunk_size;//=100;
	static int offset=0;
	static byte[] dimension = new byte[dims];
	static short[][] map;// = new short[tuples][dims];
	static LinkedList<Integer>[][] distinctValuesPerAssign;// = (LinkedList<Short>[][]) new LinkedList[chunk_size][dims];
	static int[][][] MinMaxPerAssign;
	static int[][] final_assignment;// = new int[tuples][k];
	static LinkedList chunk_sizes = new LinkedList();

	//static long timeToCalcCovers;

	//*******************************************//
	//METHODS THAT PERFORM ARRAY-PROCESSING TASKS//
	//*******************************************//

	public static double findLargest		//Finds the largest element in a positive array.
	(double[][] array)
	//works for arrays where all values are >= 0.
	{
		double largest = 0;
		for (int i=0; i<array.length; i++){
			for (int j=0; j<array[i].length; j++){
				if (array[i][j] > largest)
				{
					largest = array[i][j];
				}
			}
		}

		return largest;
	}

	public static double[][] copyOf	(double[][] original){
		double[][] copy = new double[original.length][original[0].length];
		for (int i=0; i<original.length; i++){
			//Need to do it this way, otherwise it copies only memory location
			System.arraycopy(original[i], 0, copy[i], 0, original[i].length);
		}

		return copy;
	}

	public static void transpose(double[][] cost){
		for (int i=0; i<cost.length; i++)		//Generate cost by subtracting.
		{
			for (int j=0; j<cost[i].length; j++)
			{
				cost [i][j] = (BIG - cost [i][j]);
			}
		}
	}

	
	//*******************************************//
	//METHODS to pre-process the data//
	//*******************************************//
	public static void dimension_sort() { // sort the dimensions according to their effect to GCP
		for (int i = 0;i < dims;i++) {
			dimension[i] = (byte) i;
		}
		for (int i = 0;i < dims;i++) { // small domain, use bubble sort
			for (int j = 0;j < dims-i-1;j++) { // smaller range, put it the front
				if (cardinalities[dimension[j]] > cardinalities[dimension[j+1]]) {
					byte temp = dimension[j];
					dimension[j] = dimension[j+1];
					dimension[j+1] = temp;
				}
			}
		}
	}

	// return true if data x > data y
	public static boolean compare_data(short[] x, short[] y) {
		for (int i = 0;i < dims;i++) {
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
		/*if (map[x][dimension[0]]>map[y][dimension[0]])
			System.out.println("oops");*/
		return false;
	}

	// sort data by lexicographical order


	public static void quickSort(int left, int right) {    
		int i = left, j = right;
		int ref = left + (right-left)/2;
		short[] pivot = map[ref];
		short[] temp2 = new short[dims];
		while (i <= j) {        
			while (compare_data(pivot, map[i]))
				i++;
			while (compare_data(map[j], pivot))
				j--;
			if (i <= j) {
				temp2=map[i];
				map[i]=map[j];
				map[j]=temp2;

				i++;
				j--;
			}
		};    
		// recursion
		if (left < j)
			quickSort(left, j);
		if (i < right) {
			quickSort(i, right);
		}
	}
	private static void simple_partition() {
		for (int i=0; i<tuples/chunk_size; i++)
			chunk_sizes.add(chunk_size);

	}


	static void partition(int s, int e, int dim_ref) {            
		int max = tuples-1;
		if (cardinalities[dimension[dim_ref]] + 1< tuples) {
			max = cardinalities[dimension[dim_ref]] + 1;
		}
		int[] group_size = new int[max];
		boolean[] combined_group= new boolean[max];
		// group the data

		combined_group[0] = false; // indicate whether the group is pure or not
		// if the group is not pure (it is merged with some groups), we do not recursively partition the group. We generalize the group

		int group_count = 0;
		group_size[group_count] = 1;
		for (int i = s+1;i <= e;i++) {     
			if (map[i][dimension[dim_ref]] == map[i-1][dimension[dim_ref]]) {
				group_size[group_count] ++;
			} else {
				group_count++;
				combined_group[group_count] = false;
				group_size[group_count] = 1; 
			}
		}
		group_count++;    
		// note that we just group the tuples by attribute, count may < k

		// see if each group is ok (count >= k?) and merge to neighboring group if not good
		int remain = e-s+1; // keep tract of the number of tuples remained
		for (int i = 0;i < group_count - 1;i++) {
			if (remain - group_size[i] < k) { // we have to group the remaining groups
				if (remain < 2 * k) { // merge it to the previous group
					for (int j = i;j < group_count-1;j++) {
						group_size[j] = -1; // indicate the group is removed, lazy deletion
					}
					group_size[group_count - 1] = remain;
					combined_group[group_count-1] = true;
					break;
				} else {
					for (int j = i+1;j < group_count-1;j++) {
						group_size[j] = -1;
					}
					group_size[group_count - 1] = k;
					combined_group[group_count-1] = true;
					group_size[i] = remain - k;
					break;            
				}
			} else if (group_size[i] < k) { // the group is not good
				if (i != group_count - 1 && group_size[i+1] < k) { // find somebody to merge
					group_size[i+1] += group_size[i];
					combined_group[i+1] = true;
					group_size[i] = -1;
				} else if (i == 0 || group_size[i-1] == -1) {
					if (group_size[i+1] + group_size[i] >= 2 * k) {                   
						group_size[i+1] -= k - group_size[i];
						group_size[i] = k;
						combined_group[i] = true;
					} else {
						group_size[i+1] += group_size[i];
						combined_group[i+1] = true;
						group_size[i] = -1;
					}
				} else {
					if (combined_group[i-1] || group_size[i-1] + group_size[i] < 2 * k) {
						group_size[i] += group_size[i-1];
						combined_group[i] = true;
						remain += group_size[i-1];
						group_size[i-1] = -1;    
					} else {
						remain += k - group_size[i];
						group_size[i-1] -= k - group_size[i];
						group_size[i] = k;
						combined_group[i] = true;
					}
				}
			}
			if (group_size[i] != -1) { // skip those that are already merged
				remain -= group_size[i];
			}
		} 
		int start = s;
		for (int i = 0;i < group_count;i++) {        
			if (group_size[i] != -1) {
				if (combined_group[i] || group_count == 1 || dim_ref == dims - 1) {
					chunk_sizes.addLast(group_size[i]);
					start += group_size[i];

				} else {
					partition(start, start+group_size[i]-1, dim_ref + 1);                
					start += group_size[i];
				}
			}
		}
	}

	static void partition2(int s, int e, int dim_ref) {            
		int max = cardinalities[dimension[dim_ref]] + 1;

		int[] group_size = new int[max];
		boolean[] combined_group= new boolean[max];
		// group the data

		combined_group[0] = false; // indicate whether the group is pure or not
		// if the group is not pure (it is merged with some groups), we do not recursively partition the group. We generalize the group

		int group_count = 0;
		group_size[group_count] = 1;
		for (int i = s+1;i <= e;i++) {     
			if (map[i][dimension[dim_ref]] == map[i-1][dimension[dim_ref]]) {
				group_size[group_count] ++;
			} else {
				group_count++;
				combined_group[group_count] = false;
				group_size[group_count] = 1; 
			}
		}
		group_count++;    
		// note that we just group the tuples by attribute, count may < k

		// see if each group is ok (count >= k?) and merge to neighboring group if not good
		int remain = e-s+1; // keep tract of the number of tuples remained
		for (int i = 0;i < group_count - 1;i++) {
			if (remain - group_size[i] < partition_size) { // we have to group the remaining groups
				if (remain < 2 * partition_size) { // merge it to the previous group
					for (int j = i;j < group_count-1;j++) {
						group_size[j] = -1; // indicate the group is removed, lazy deletion
					}
					group_size[group_count - 1] = remain;
					combined_group[group_count-1] = true;
					break;
				} else {
					for (int j = i+1;j < group_count-1;j++) {
						group_size[j] = -1;
					}
					group_size[group_count - 1] = remain - group_size[i];
					combined_group[group_count-1] = true;
					//group_size[i] = remain - k;
					break;            
				}
			} else if (group_size[i] < partition_size) { // the group is not good
				if (i != group_count - 1 && group_size[i+1] < partition_size) { // find somebody to merge
					group_size[i+1] += group_size[i];
					combined_group[i+1] = true;
					group_size[i] = -1;
				} else if (i == 0 || group_size[i-1] == -1) {
					if (group_size[i+1] + group_size[i] >= 2 * partition_size) {                   
						/*group_size[i+1] -= k - group_size[i];
						group_size[i] = k;
						combined_group[i] = true;*/
					} else {
						group_size[i+1] += group_size[i];
						combined_group[i+1] = true;
						group_size[i] = -1;
					}
				} else {
					if (combined_group[i-1] || group_size[i-1] + group_size[i] < 2 * k) {
						group_size[i] += group_size[i-1];
						combined_group[i] = true;
						remain += group_size[i-1];
						group_size[i-1] = -1;    
					} else {
						/*remain += k - group_size[i];
						group_size[i-1] -= k - group_size[i];
						group_size[i] = k;
						combined_group[i] = true;*/
					}
				}
			}
			if (group_size[i] != -1) { // skip those that are already merged
				remain -= group_size[i];
			}
		} 
		int start = s;
		for (int i = 0;i < group_count;i++) {        
			if (group_size[i] != -1) {
				if (combined_group[i] || group_count == 1 || dim_ref == dims - 1) {
					chunk_sizes.addLast(group_size[i]);
					start += group_size[i];

				} else {
					partition2(start, start+group_size[i]-1, dim_ref + 1);                
					start += group_size[i];
				}
			}
		}
	}


	//***********//
	//MAIN METHOD//
	//***********//

	public static void main(String[] args) 	{
		//timeToCalcCovers = 0;
		TopCoderAlgo algo = new TopCoderAlgo();
		String inputFile = args[0];
		partition_size = Integer.parseInt(args[1]);
		k = Integer.parseInt(args[2]);
		tuples = Integer.parseInt(args[3]);
		chunk_size = Integer.parseInt(args[1]);

		map = new short[tuples][dims];
		final_assignment = new int[tuples][k];

		CensusParser tp = new CensusParser(inputFile, dims);
		System.out.println(inputFile+" "+partition_size +" "+ k + "  "+ tuples);
		try {
			int i=0;
			while (tp.hasNext()){
				map[i++]=tp.nextTuple2();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		//Below enter "max" or "min" to find maximum sum or minimum sum assignment.
		String sumType = "min";		

		//sort tuples
		dimension_sort();
		quickSort(0, tuples-1);
		//partition(0, tuples-1, 0);
		//partition2(0,tuples-1, 0);
		simple_partition();
		long startTime = System.currentTimeMillis();
		long time_of_hungarian = 0;
		long start_of_hungarian = 0;
		int index = 0;
		chunk_size = (Integer) chunk_sizes.get(index);
		double distortion = 0.0;
		while (offset<tuples){
			if (MIXED){
				MinMaxPerAssign = new int[chunk_size][dims][2];
				distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];
			}else
				if (RANGE)
					MinMaxPerAssign = new int[chunk_size][dims][2];
				else
					distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];
			//Hard-coded example.
			double[][] array = computeCostMatrix();

			int[] assignment = new int[array.length];

			int times = 0;

			while (++times<k){
				//hgAlgorithm(array, sumType, assignment);	//Call Hungarian algorithm.
				start_of_hungarian = System.currentTimeMillis();
				algo.hungarian(array, assignment);
				time_of_hungarian+=(System.currentTimeMillis() - start_of_hungarian);

                System.out.println("time "+times);
				for (int i=0; i<assignment.length; i++){
					//final_assignment[assignment[i][0]+offset][times]= assignment[i][1]+offset;
					//findSet(assignment[i][0], map[assignment[i][1]+offset]);
					final_assignment[offset+i][times]= assignment[i]+offset;
					if (MIXED){
						findSet_mixed(i, map[assignment[i]+offset]);
					}else{
						if (RANGE)
							findSet_numerical(i, map[assignment[i]+offset]);
						else
							findSet(i, map[assignment[i]+offset]);
					}
				}
				if (times!=k-1)
					recomputeCostMatrix(array, times);
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

            //**** BEGIN XUE MINGQIANG **** //
            // this call returns a random assignment generated from the k-regular matching graph
            int [] rand_A = Randomization.run(final_assignment, offset, assignment.length, k); 
            //**** END XUE MINGQIANG **** //

			
			/*if (OPTIMIZATION){
				double previousDistortion = distortion;
				distortion = 0.0;
				while (distortion<previousDistortion){
					distortion = 0.0;
					array = OptimizationCostMatrix();
				}
			}*/
			

			offset+=chunk_size;
			if (index<chunk_sizes.size()-1)
				chunk_size = (Integer) chunk_sizes.get(++index);
		}

		long endTime = System.currentTimeMillis();

		System.out.println("The winning assignment after "+index+" runs (" + sumType + " sum) is:\n");	


		for (int i=0; i<final_assignment.length; i++){
			for (int j=0; j<k; j++){
				System.out.print(final_assignment[i][j]+" ");
			}
			//distortion+= NCP(distinctValuesPerAssign[i]);
			//distortion+= array[i][assignment[i][1]];
			System.out.println();

		}

		


		System.out.println("Time: "+(endTime - startTime)+"  "+time_of_hungarian+"\n Distortion "+ (double)(distortion/(dims*tuples)));
		//System.out.println(timeToCalcCovers);
	}


	private static double[][] computeCostMatrix() {
		int totalTrunc = 0;
		double[][] cost = new double[chunk_size][chunk_size];
		for (int i=0; i<chunk_size; i++){
			for (int j=i; j<chunk_size; j++){

				if (i==j)
					cost[i][j]=BIG;
				else{
					double c;
					if (MIXED){
						c=NCP_mixed(map[offset+i], map[offset+j]);
					}else
						if (RANGE)
							c=NCP_numerical(map[offset+i], map[offset+j]);
						else
							c=NCP(map[offset+i], map[offset+j]);
					if (c<=maxCost){
						cost[i][j]=c;
						cost[j][i]=c;
					}else{
						totalTrunc++;
						cost[i][j]=BIG;
						cost[j][i]=BIG;
					}
				}
			}
			//System.out.println(totalTrunc);
			final_assignment[offset+i][0]= offset+i;
			if (MIXED){
				for (int l=0; l<dims; l++){
					if (l==0||l==2){
						MinMaxPerAssign[i][l][0]=map[offset+i][l];
						MinMaxPerAssign[i][l][1]=map[offset+i][l];
					}else{
						LinkedList<Integer> list = new LinkedList<Integer>();
						list.add((int) map[offset+i][l]);
						distinctValuesPerAssign[i][l] = list;
					}

				}
			}else
				if (RANGE){
					for (int l=0; l<dims; l++){
						MinMaxPerAssign[i][l][0]=map[offset+i][l];
						MinMaxPerAssign[i][l][1]=map[offset+i][l];
					}
				}else{
					for (int l=0; l<dims; l++){
						LinkedList<Integer> list = new LinkedList<Integer>();
						list.add((int) map[offset+i][l]);
						distinctValuesPerAssign[i][l] = list;
					}
				}
		}

		return cost;
	}

	private static void recomputeCostMatrix(double[][] array, int times) {
		for (int i=0; i<chunk_size; i++){
			int[] tuples_per_assignment = final_assignment[offset+i];
			for (int l=0; l<=times; l++)
				array[i][tuples_per_assignment[l]-offset]=BIG;

			for (int j=0; j<chunk_size; j++){
				if (array[i][j]==BIG)
					continue;
				else {
					if (MIXED){
						array[i][j]=NCP_mixed(map[offset+j], MinMaxPerAssign[i],distinctValuesPerAssign[i]);
					}else
						if (RANGE)
							array[i][j]=NCP_numerical(map[offset+j], MinMaxPerAssign[i]);
						else
							array[i][j]=NCP(map[offset+j], (distinctValuesPerAssign[i]));
				}

			}
		}

	}
	///////////////////Mixed Representation//////////////	
	private static double NCP_mixed(short[] tuple1, short[] tuple2){
		double score=0.0;
		for (int i=0; i<dims; i++){
			if (i==0 || i==2)
				score+=(double)Math.abs(tuple1[i]-tuple2[i])/(double)(cardinalities[i]-1);
			else{
				if (tuple1[i]==tuple2[i])
					score+=0;
				else 
					score+=(double)(1)/(double)(cardinalities[i]-1);
			}	
		}
		return score;
	}

	private static double NCP_mixed(short[] tuple, int[][] MinMaxPerDim, LinkedList<Integer>[] distinctValuesPerDim ){
		double score=0.0;
		int min;
		int max;
		for (int i=0; i<dims; i++){
			if (i==0 || i==2){
				int[] distinctValues = MinMaxPerDim[i];
				min = distinctValues[0];
				max = distinctValues[1];
				if (tuple[i] < min)
					score+=(double)(max-tuple[i])/(double)(cardinalities[i]-1);
				else if (tuple[i]>max)
					score+=(double)(tuple[i]-min)/(double)(cardinalities[i]-1);
				else
					score+=(double)(max-min)/(double)(cardinalities[i]-1);
			}else{
				LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
				if (!distinctValues.contains((int)tuple[i]))
					score+=(double)(distinctValues.size())/(double)(cardinalities[i]-1);
				else 
					score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
			}
		}
		return score;
	}
	private static double NCP_mixed(int[][] MinMaxPerDim, LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;
		for (int i=0; i<dims; i++){
			if (i==0 || i==2){
				int[] distinctValues = MinMaxPerDim[i];
				score+=(double)(distinctValues[1]-distinctValues[0])/(double)(cardinalities[i]-1);
			}else{
				LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
				score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
			}
		}
		return score;
	}

	private static void findSet_mixed(int assign_number, short[] newTuple){
		LinkedList<Integer>[] dValues = distinctValuesPerAssign[assign_number];
		int[][] MinMaxPerDim = MinMaxPerAssign[assign_number];
		for (int i=0; i<dims; i++){
			if (i==0||i==2){
				int[] distinctValues = MinMaxPerDim[i];
				if (newTuple[i] < distinctValues[0]){
					distinctValues[0]=(int) newTuple[i];
				}
				else if (newTuple[i] > distinctValues[1]){
					distinctValues[1]=(int) newTuple[i];
				}
			}else{
				LinkedList distValuesPerDim = dValues[i];
				if (!distValuesPerDim.contains((int)newTuple[i]))
					distValuesPerDim.add((int)newTuple[i]);
			}

		}
		return;
	}


	///////////////////Range Representation//////////////			
	private static double NCP_numerical(short[] tuple1, short[] tuple2){
		double score=0.0;
		for (int i=0; i<dims; i++){
			score+=(double)Math.abs(tuple1[i]-tuple2[i])/(double)(cardinalities[i]-1);

		}
		return score;
	}
	private static double NCP_numerical(short[] tuple, int[][] MinMaxPerDim){
		double score=0.0;
		int min;
		int max;
		for (int i=0; i<dims; i++){
			int[] distinctValues = MinMaxPerDim[i];
			min = distinctValues[0];
			max = distinctValues[1];
			if (tuple[i] < min)
				score+=(double)(max-tuple[i])/(double)(cardinalities[i]-1);
			else if (tuple[i]>max)
				score+=(double)(tuple[i]-min)/(double)(cardinalities[i]-1);
			else
				score+=(double)(max-min)/(double)(cardinalities[i]-1);
		}
		return score;
	}
	private static double NCP_numerical(int[][] MinMaxPerDim){
		double score=0.0;
		for (int i=0; i<dims; i++){
			int[] distinctValues = MinMaxPerDim[i];
			score+=(double)(distinctValues[1]-distinctValues[0])/(double)(cardinalities[i]-1);
		}
		return score;
	}
	private static void findSet_numerical(int assign_number, short[] newTuple){
		int[][] MinMaxPerDim = MinMaxPerAssign[assign_number];
		for (int i=0; i<MinMaxPerDim.length; i++){
			int[] distinctValues = MinMaxPerDim[i];
			if (newTuple[i] < distinctValues[0]){
				distinctValues[0]=(int) newTuple[i];
			}
			else if (newTuple[i] > distinctValues[1]){
				distinctValues[1]=(int) newTuple[i];
			}

		}
		return;
	}

	////////////////Set Representation///////////////
	private static double NCP(short[] tuple1, short[] tuple2){
		double score=0.0;
		for (int i=0; i<dims; i++){
			if (tuple1[i]==tuple2[i])
				score+=0;
			else 
				score+=(double)(1)/(double)(cardinalities[i]-1);
		}
		return score;
	}

	private static double NCP(short[] tuple, LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;
		for (int i=0; i<dims; i++){
			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			if (!distinctValues.contains((int)tuple[i]))
				score+=(double)(distinctValues.size())/(double)(cardinalities[i]-1);
			else 
				score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);
		}
		return score;
	}
	private static double NCP(LinkedList<Integer>[] distinctValuesPerDim){
		double score=0.0;
		for (int i=0; i<dims; i++){
			LinkedList<Integer> distinctValues = distinctValuesPerDim[i];
			score+=(double)(distinctValues.size()-1)/(double)(cardinalities[i]-1);

		}
		return score;
	}

	private static void findSet(int assign_number, short[] newTuple){
		LinkedList<Integer>[] distinctValues = distinctValuesPerAssign[assign_number];
		for (int i=0; i<distinctValues.length; i++){
			LinkedList distValuesPerDim = distinctValues[i];
			if (!distValuesPerDim.contains((int)newTuple[i]))
				distValuesPerDim.add((int)newTuple[i]);
		}
		return;
	}
}
