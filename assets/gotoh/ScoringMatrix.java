import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


public class ScoringMatrix {

	private char[] indexToChar;
	private int[] charToIndex;
	private int[][] svalues;
	private double[][] values;
	private int mult=0;
	
	
	public int[] encodeSequence(String sequence) {
		int[] re = new int[sequence.length()];
		for (int i=0; i<sequence.length(); i++)
			re[i] = charToIndex[sequence.charAt(i)];
		return re;
	}
	
	public char[] decodeSequence(int[] sequence) {
		char[] re = new char[sequence.length];
		for (int i=0; i<sequence.length; i++)
			re[i] = indexToChar[sequence[i]];
		return re;
	}
	
	public char[] getIndexToChar() {
		return indexToChar;
	}
	
	public int[][] getScaledMatrix() {
		return svalues;
	}
	public double[][] getMatrix() {
		return values;
	}
	
	public int[] getCharToIndex() {
		return charToIndex;
	}
	
	public int getMult() {
		return mult;
	}
	
	private void setChars(String chars){
		indexToChar = chars.toCharArray();
		charToIndex = new int[256];
		for (int i=0; i<chars.length(); i++) 
			charToIndex[chars.charAt(i)] = i;
	}
	
	private void setMatrix(double[][] matrix, int sig) {
		mult = (int)Math.pow(10, sig);
		svalues = new int[matrix.length][matrix.length];
		values = new double[matrix.length][matrix.length];
		for (int i=0; i<matrix.length; i++)
			for (int j=0; j<matrix.length; j++) {
				svalues[i][j] = (int)Math.round(matrix[i][j]*mult);
				values[i][j] = matrix[i][j];
			}
	}
	
	public static ScoringMatrix fromFile(File file) throws IOException {
		ScoringMatrix re = new ScoringMatrix();
		String chars = null;
		double[][] matrix = null;
		int nMatrix = 0;
		int sig = 0;
		boolean full = false;
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while ((line=br.readLine())!=null) {
			if (line.startsWith("ROWINDEX") || line.startsWith("COLINDEX")) {
				if (chars==null) {
					chars = line.split("\\s+")[1];
					re.setChars(chars);
					matrix = new double[chars.length()][chars.length()];
				} else {
					if (!chars.equals(line.split("\\s+")[1]))
						throw new IOException("ROWINDEX and COLINDEX must be equal!");
				}
			}
			else if (line.startsWith("MATRIX")) {
				if (chars==null)
					throw new IOException("ROWINDEX and COLINDEX be specified before MATRIX!");
				
				String[] numbers = line.split("\\s+");
				full |= nMatrix==0 && numbers.length-1>1;
				for (int i=0; i<numbers.length-1; i++) {
					
					sig = Math.max(sig, numbers[i+1].indexOf('.')>=0 ? numbers[i+1].length()-numbers[i+1].indexOf('.')-1: 0);
					
					matrix[nMatrix][i] = Double.parseDouble(numbers[i+1]);
//					System.out.println(chars.charAt(nMatrix)+" "+chars.charAt(i)+" -> "+matrix[nMatrix][i]);
					if (!full)
						matrix[i][nMatrix] = matrix[nMatrix][i]; 
				}
				
				nMatrix++;
			}
		}
		br.close();
		re.setMatrix(matrix,sig);
		return re;		
	}
	
	
}
