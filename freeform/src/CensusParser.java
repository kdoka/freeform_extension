import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;


public class CensusParser {
	private DataInputStream dis;
	int dims;
	
	public CensusParser(String tuplesFile ) {
		//this.dims = dims;
		try {
			dis = new DataInputStream(new BufferedInputStream(new FileInputStream(tuplesFile)));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	public CensusParser(String tuplesFile, int dims ) {
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
	
	public short[] nextTuple2() throws IOException{
		short[] tuple =  new short[dims];
		String line = dis.readLine();
		String[] splitLine = line.split(", ");
		for (int i=0; i<dims;i++){
			tuple[i]=(short) Integer.parseInt(splitLine[i].trim());
		}
		
		return tuple;
	}
}
