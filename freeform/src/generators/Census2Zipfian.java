package generators;

import java.io.IOException;
import java.util.*;
import java.io.FileWriter;

public class Census2Zipfian {

	static int SA = 7;// 0 - 7.
	static int dims = 8; //3
	static int tuples;// = 10000;
	static int[][] map;// = new int[tuples][dims];
	static int minValue = 100;
	static int maxValue = 0;
	static double theta;

	
	public static void main(String[] args) 	{

		if (args.length!=3){
			System.out.println("\nUsage:   java Census2Zipfian inFile #tuples theta");
			System.out.println("Changes the last attribute to a synthetic zipfian distribution.");
			System.out.println("min and max values are the same as in the original file.");
//			System.out.println("\t SA: index of sensitive attribute [0-7].");
			return;
		}

		String inputFile = args[0];
		tuples = Integer.parseInt(args[1]);
		theta = Double.parseDouble(args[2]);
//		SA = Integer.parseInt(args[3]);

		map = new int[tuples][dims];
		
		try {
			CensusParserInt tp = new CensusParserInt(inputFile, dims);
			int i=0;
			while (tp.hasNext()){
				map[i]=tp.nextTuple2();
				if(map[i][SA] < minValue){
					minValue = map[i][SA];
				}
				if(map[i][SA] > maxValue){
					maxValue = map[i][SA];
				}
				i++;
			}
		} catch (IOException ioe) {
			ioe.printStackTrace();
		} catch (Exception e){
			e.printStackTrace();
		}
		System.out.println("min="+minValue+" max="+maxValue);
		int range = (maxValue - minValue) + 1; //integers.
		Stats gen = new Stats(range, theta);

//		for (int i=0; i<map.length; i++){
//			map[i][SA] = gen.zipfGen();
//			map[i][SA] += minValue;
//		}
		
		//Save Results:
		FileWriter fw = null;
		try{
			String outfile = "./zipfDataU"+theta+"_"+tuples+".txt";
			fw = new FileWriter(outfile,true); //true == append
			for (int i=0; i<map.length; i++){
				for (int j=0; j<dims-1; j++){
					fw.write(map[i][j]+", ");
				}
				int newVal = minValue + gen.zipfGen();
				fw.write(newVal+"\n");
			}
		}catch(IOException ioe){
			System.err.println("IOException: " + ioe.getMessage());
		}finally{
			try{
				if(fw != null) fw.close();
			}catch(Exception e){
				//ignore.
			}
		}

	}

}
