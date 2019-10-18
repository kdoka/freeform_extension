import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

public class SortGreedy {

	private static final double BIG = Double.MAX_VALUE;
	private static final boolean RANGE = false;
	private static final boolean MIXED = true;
	private static double maxCost = BIG;//0.5436267458828434;
	static byte[] cardinalities = {79, 2, 17, 6, 9, 10, 83, 51};
	//static byte[] cardinalities ={41, 41};
	static int partition_size;// = 1000;
	static int k;// = 10;
	static int dims;
	static int tuples;// = 10000;
	static int chunk_size;//=100;
	static int offset=0;
	static byte[] dimension;
	//static short[][] map = {{59,25}, {57, 27}, {39, 47}, {28, 41}, {41, 20}, {37, 59}, {40, 35}, {53, 34}};

	static short[][] map;// = new short[tuples][dims];
	//static short[][] map = {{43, 37}, {26, 29}, {42, 31}, {32, 48}, {48, 60}, {44, 27}, {40, 41}, {39, 25}};
	//static short[][] map = {{44, 40}, {34, 32}, {50, 59}, {41, 26}, {57, 36}, {33, 50}, {29, 35}, {26, 45}};
	//static short[][] map = {{56, 57}, {33, 60}, {46, 20}, {20, 21}, {49, 53}, {26, 44}, {55, 42}, {37, 56}};
	//static short[][] map = {{52, 38}, {38, 33}, {23, 51}, {31, 30}, {20, 40}, {47, 55}, {34, 39}, {50, 35}};
	//static short[][] map = {{50, 41}, {59, 22}, {60, 22}, {45, 31}, {50, 20}, {20, 54}, {53, 46}, {55, 31}};

	static LinkedList<Integer>[][] distinctValuesPerAssign;// = (LinkedList<Short>[][]) new LinkedList[chunk_size][dims];
	static int[][][] MinMaxPerAssign;
	static int[][] final_assignment;// = new int[tuples][k];
	static LinkedList chunk_sizes = new LinkedList();
	static HeapNode[] edges;
	static int edge_size;

	static long backtrace_time;
	static long cost_matrix_time = 0;;
	static double threshold;

	//*******************************************//
	//METHODS THAT PERFORM ARRAY-PROCESSING TASKS//
	//*******************************************//
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

	private static int[] greedyAssign(double[][] array, int[] assignment) {

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
				return true;
			} else if (x[dimension[i]] < y[dimension[i]]) {
				return false;
			}        
		} 
		return false;
	}

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
		backtrace_time = 0;
		//timeToCalcCovers = 0;
		//HungarianAlgorithm algo= new HungarianAlgorithm();
		String inputFile = args[0];
		partition_size = Integer.parseInt(args[1]);
		k = Integer.parseInt(args[2]);
		tuples = Integer.parseInt(args[3]);
		chunk_size = Integer.parseInt(args[1]);
		dims = Integer.parseInt(args[4]);
		dimension = new byte[dims];
		map = new short[tuples][dims];
		final_assignment = new int[tuples][k];
		edges = new HeapNode[tuples*tuples];

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
		long delta = 0;
		long t_start;
		int index = 0;
		chunk_size = (Integer) chunk_sizes.get(index);
		double distortion = 0.0;
		long edge_sort_time = 0;
		long edge_sort_start;
		while (offset<tuples){
			if (MIXED){
				MinMaxPerAssign = new int[chunk_size][dims][2];
				distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];
			}else
				if (RANGE)
					MinMaxPerAssign = new int[chunk_size][dims][2];
				else
					distinctValuesPerAssign = (LinkedList<Integer>[][]) new LinkedList[chunk_size][dims];

			
			double[][] array = computeCostMatrix();

			int[] assignment = new int[array.length];

			int times = 0;

			edge_sort_start = System.currentTimeMillis();
			qSort(0, edge_size-1);
			//Collections.sort(edges, new HeapComparator());
			edge_sort_time += System.currentTimeMillis() - edge_sort_start;


			while (++times<k){


				long backtrace_start = System.currentTimeMillis();
				greedyAssign(array, assignment);	//Call Hungarian algorithm.
				backtrace_time += (System.currentTimeMillis()-backtrace_start);
				//assignment=algo.computeAssignments(array);
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
				if (times!=k-1){
					recomputeCostMatrix2(array, times);
					edge_sort_start = System.currentTimeMillis();
					qSort(0, edge_size-1);
					//Collections.sort(edges, new HeapComparator());
					edge_sort_time += System.currentTimeMillis() - edge_sort_start;
				}
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
			t_start = System.currentTimeMillis();
			int [] rand_A = Randomization.run(final_assignment, offset, assignment.length, k); 
			delta += System.currentTimeMillis() - t_start;
			//**** END XUE MINGQIANG **** //

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

		System.out.println("Time: "+(endTime - startTime)+"\n Distortion "+ (double)(distortion/(dims*tuples)));
		System.out.println("Time without randomization: "+(endTime - startTime - delta));
		System.out.println("Time for sorting the edge list: "+edge_sort_time);
		System.out.println("Backtrace time: "+backtrace_time);
		System.out.println("Cost Matrix time: "+cost_matrix_time);
		//System.out.println(timeToCalcCovers);
	}


	private static double[][] computeCostMatrix() {
		long cost_start = System.currentTimeMillis();
		edge_size = 0;
		double[][] cost = new double[chunk_size][chunk_size];
		for (int i=0; i<chunk_size; i++){
			for (int j=i; j<chunk_size; j++){

				if (i==j){

					cost[i][j]=BIG;
					/*HeapNode hn = new HeapNode(i,j,BIG);
					edges[index++]=hn;*/
				}
				else{
					double c;
					if (MIXED){
						c=NCP_mixed(map[offset+i], map[offset+j]);
					}else
						if (RANGE)
							c=NCP_numerical(map[offset+i], map[offset+j]);
						else
							c=NCP(map[offset+i], map[offset+j]);

					cost[i][j]=c;
					cost[j][i]=c;


					HeapNode hn = new HeapNode(i,j,c);
					HeapNode hn2 = new HeapNode(j,i,c);
					edges[edge_size++]=hn;
					edges[edge_size++]=hn2;
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
		cost_matrix_time+=System.currentTimeMillis()-cost_start;
		return cost;
	}

	
	private static void recomputeCostMatrix2(double[][] array, int times) {
		long cost_start = System.currentTimeMillis();
		edge_size=0;
		for (int i=0; i<chunk_size; i++){
			int[] tuples_per_assignment = final_assignment[offset+i];
			for (int l=0; l<=times; l++){
				array[i][tuples_per_assignment[l]-offset]=BIG;
				//HeapNode hn = new HeapNode(i,tuples_per_assignment[l]-offset,BIG);
				//edges[i*chunk_size+tuples_per_assignment[l]-offset]=hn;
			}

			for (int j=0; j<chunk_size; j++){
				if (array[i][j]==BIG){
					;
					//HeapNode hn = new HeapNode(i,j,BIG);
					//edges[i*chunk_size+j]=hn;
				}
				else {
					if (MIXED){
						array[i][j]=NCP_mixed(map[offset+j], MinMaxPerAssign[i],distinctValuesPerAssign[i]);
						HeapNode hn = new HeapNode(i,j,NCP_mixed(map[offset+j], MinMaxPerAssign[i],distinctValuesPerAssign[i]));
						edges[edge_size++]=hn;
					}else
						if (RANGE){
							array[i][j]=NCP_numerical(map[offset+j], MinMaxPerAssign[i]);
							HeapNode hn = new HeapNode(i,j,NCP_numerical(map[offset+j], MinMaxPerAssign[i]));
							edges[edge_size++]=hn;
						}
						else{
							array[i][j]=NCP(map[offset+j], (distinctValuesPerAssign[i]));
							HeapNode hn = new HeapNode(i,j,NCP(map[offset+j], (distinctValuesPerAssign[i])));
							edges[edge_size++]=hn;
						}
				}

			}
		}
		cost_matrix_time+=System.currentTimeMillis()-cost_start;

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
