package generators;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;


public class CensusParserInt {
	private DataInputStream dis;
	int dims;
	
	public CensusParserInt(String tuplesFile ) {
		//this.dims = dims;
		try {
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream(tuplesFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	public CensusParserInt(String tuplesFile, int dims ) {
		this.dims = dims;
		try {
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream(tuplesFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public boolean hasNext() throws IOException{
		if (dis.available()>1)
			return true;
		dis.close();
		return false;
	}
	public String nextTuple()throws IOException{
		return dis.readLine();
		
	}
	
	public int[] nextTuple2() throws IOException{
		int[] tuple =  new int[dims];
		String line = dis.readLine();
		String[] splitLine = line.split(", ");
		for (int i=0; i<dims;i++){
			tuple[i]= Integer.parseInt(splitLine[i]);
		}
		
		return tuple;
	}
}
