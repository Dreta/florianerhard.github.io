import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;


/**
 * This is quite ugly code, I know! It is the product of 7 year old optimized code with deoptimizations...
**/
public class Gotoh {

	private static final int INF = Integer.MIN_VALUE/2;

	private int sgapOpen;
	private int sgapExtend;
	private double gapOpen;
	private double gapExtend;

	private HashMap<String,Integer> sscoringMatrix1;
	private HashMap<Character,HashMap<Character,Integer>> sscoringMatrix2;
	private int[][] sscoringMatrix3;
	private int[][] sscoringMatrix4;
	
	private HashMap<String,Double> scoringMatrix1;
	private HashMap<Character,HashMap<Character,Double>> scoringMatrix2;
	private double[][] scoringMatrix3;
	private double[][] scoringMatrix4;
	
	private float mult;


	public Gotoh(int maxLength, ScoringMatrix m, int gapOpen, int gapExtend) {
		this.sscoringMatrix4 = m.getScaledMatrix();
		this.scoringMatrix4 = m.getMatrix();
		
		this.scoringMatrix3 = new double[256][256];
		this.sscoringMatrix3 = new int[256][256];
		this.scoringMatrix2 = new HashMap<Character, HashMap<Character,Double>>();
		this.sscoringMatrix2 = new HashMap<Character, HashMap<Character,Integer>>();
		this.scoringMatrix1 = new HashMap<String, Double>();
		this.sscoringMatrix1 = new HashMap<String, Integer>();
		
		char[] a = m.getIndexToChar();
		for (int i1=0; i1<a.length; i1++) {
			for (int i2=0; i2<a.length; i2++) {
				this.scoringMatrix3[a[i1]][a[i2]] = this.scoringMatrix4[i1][i2];
				this.sscoringMatrix3[a[i1]][a[i2]] = this.sscoringMatrix4[i1][i2];
				this.scoringMatrix2.computeIfAbsent(a[i1],x->new HashMap<>()).put(a[i2],this.scoringMatrix4[i1][i2]);
				this.sscoringMatrix2.computeIfAbsent(a[i1],x->new HashMap<>()).put(a[i2],this.sscoringMatrix4[i1][i2]);
				this.scoringMatrix1.put(new String(new char[] {a[i1],a[i2]}), this.scoringMatrix4[i1][i2]);
				this.sscoringMatrix1.put(new String(new char[] {a[i1],a[i2]}), this.sscoringMatrix4[i1][i2]);
			}	
		}
		
		
		this.sgapOpen = gapOpen*m.getMult();
		this.sgapExtend = gapExtend*m.getMult();
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.mult = m.getMult();
		init(maxLength,maxLength);
		initScaled(maxLength,maxLength);
	}
	
	private void init(int n, int m) {
		M = new double[n+2][m+2];
		I = new double[n+2][m+2];
		D = new double[n+2][m+2];

		I[0][0] = 0;
		D[0][0] = 0;
		M[0][0] = 0;

		double gapFirst = this.gapOpen+this.gapExtend;

		boolean initZero = false;

		M[0][1] = initZero?0:gapFirst;
		I[0][1] = INF;

		D[1][0] = INF;
		M[1][0] = initZero?0:gapFirst;

		for (int j=2; j<=m; j++) {
			M[0][j] = initZero?0:M[0][j-1]+this.gapExtend;
			I[0][j] = INF;
		}

		for (int i=2; i<=n; i++) {
			D[i][0] = INF;
			M[i][0] = initZero?0:M[i-1][0]+this.gapExtend;
		}
	}

	private void initScaled(int n, int m) {
		sM = new int[n+2][m+2];
		sI = new int[n+2][m+2];
		sD = new int[n+2][m+2];

		sI[0][0] = 0;
		sD[0][0] = 0;
		sM[0][0] = 0;

		int gapFirst = this.sgapOpen+this.sgapExtend;

		boolean initZero = false;

		sM[0][1] = initZero?0:gapFirst;
		sI[0][1] = INF;

		sD[1][0] = INF;
		sM[1][0] = initZero?0:gapFirst;

		for (int j=2; j<=m; j++) {
			sM[0][j] = initZero?0:sM[0][j-1]+this.sgapExtend;
			sI[0][j] = INF;
		}

		for (int i=2; i<=n; i++) {
			sD[i][0] = INF;
			sM[i][0] = initZero?0:sM[i-1][0]+this.sgapExtend;
		}
	}


	private int[][] sM;
	private int[][] sI;
	private int[][] sD;

	private double[][] M;
	private double[][] I;
	private double[][] D;

	public double align(String s1, String s2, int[] s1e, int[] s2e, int mode, boolean init) {
		final int n =s1e.length;
		final int m = s2e.length;

		final double gapExtend = this.gapExtend;
		final double gapFirst = this.gapOpen+this.gapExtend;
		if (init) 
			init(s1e.length, s2e.length);

		switch (mode) {
		case 1:
			for (int i=1; i<=n; i++) {
				for (int j=1; j<=m; j++) {
					I[i][j] = Math.max(I[i-1][j]+gapExtend, M[i-1][j]+gapFirst);
					D[i][j] = Math.max(D[i][j-1]+gapExtend, M[i][j-1]+gapFirst);
					String s = String.valueOf(new char[] {s1.charAt(i-1),s2.charAt(j-1)});
					M[i][j] = Math.max(M[i-1][j-1]+scoringMatrix1.get(s), Math.max(I[i][j], D[i][j]));
				}
			}
			
			break;
		case 2:
			for (int i=1; i<=n; i++) {
				HashMap<Character, Double> s = scoringMatrix2.get(s1.charAt(i-1));
				for (int j=1; j<=m; j++) {
					I[i][j] = Math.max(I[i-1][j]+gapExtend, M[i-1][j]+gapFirst);
					D[i][j] = Math.max(D[i][j-1]+gapExtend, M[i][j-1]+gapFirst);
					M[i][j] = Math.max(M[i-1][j-1]+s.get(s2.charAt(j-1)), Math.max(I[i][j], D[i][j]));
				}
			}
			
			break;
		case 3:
			for (int i=1; i<=n; i++) {
				double[] s = scoringMatrix3[s1.charAt(i-1)];
				for (int j=1; j<=m; j++) {
					I[i][j] = Math.max(I[i-1][j]+gapExtend, M[i-1][j]+gapFirst);
					D[i][j] = Math.max(D[i][j-1]+gapExtend, M[i][j-1]+gapFirst);
					M[i][j] = Math.max(M[i-1][j-1]+s[s2.charAt(j-1)], Math.max(I[i][j], D[i][j]));
				}
			}
			
			break;
		case 4:
			for (int i=1; i<=n; i++) {
				double[] s = scoringMatrix4[s1e[i-1]];
				for (int j=1; j<=m; j++) {
					I[i][j] = Math.max(I[i-1][j]+gapExtend, M[i-1][j]+gapFirst);
					D[i][j] = Math.max(D[i][j-1]+gapExtend, M[i][j-1]+gapFirst);
					M[i][j] = Math.max(M[i-1][j-1]+s[s2e[j-1]], Math.max(I[i][j], D[i][j]));
				}
			}
			break;
		}
		return M[n][m];
	}
	
	public float alignScaled(String s1, String s2, int[] s1e, int[] s2e, int mode, boolean init) {
		final int n =s1e.length;
		final int m = s2e.length;

		final int gapExtend = this.sgapExtend;
		final int gapFirst = this.sgapOpen+this.sgapExtend;
		if (init) 
			initScaled(s1e.length, s2e.length);
		
		switch (mode) {
		case 1:
			for (int i=1; i<=n; i++) {
				for (int j=1; j<=m; j++) {
					sI[i][j] = Math.max(sI[i-1][j]+gapExtend, sM[i-1][j]+gapFirst);
					sD[i][j] = Math.max(sD[i][j-1]+gapExtend, sM[i][j-1]+gapFirst);
					String s = String.valueOf(new char[] {s1.charAt(i-1),s2.charAt(j-1)});
					sM[i][j] = Math.max(sM[i-1][j-1]+sscoringMatrix1.get(s), Math.max(sI[i][j], sD[i][j]));
				}
			}
			
			break;
		case 2:
			for (int i=1; i<=n; i++) {
				HashMap<Character, Integer> s = sscoringMatrix2.get(s1.charAt(i-1));
				for (int j=1; j<=m; j++) {
					sI[i][j] = Math.max(sI[i-1][j]+gapExtend, sM[i-1][j]+gapFirst);
					sD[i][j] = Math.max(sD[i][j-1]+gapExtend, sM[i][j-1]+gapFirst);
					sM[i][j] = Math.max(sM[i-1][j-1]+s.get(s2.charAt(j-1)), Math.max(sI[i][j], sD[i][j]));
				}
			}
			
			break;
		case 3:
			for (int i=1; i<=n; i++) {
				int[] s = sscoringMatrix3[s1.charAt(i-1)];
				for (int j=1; j<=m; j++) {
					sI[i][j] = Math.max(sI[i-1][j]+gapExtend, sM[i-1][j]+gapFirst);
					sD[i][j] = Math.max(sD[i][j-1]+gapExtend, sM[i][j-1]+gapFirst);
					sM[i][j] = Math.max(sM[i-1][j-1]+s[s2.charAt(j-1)], Math.max(sI[i][j], sD[i][j]));
				}
			}
			
			break;
		case 4:
			for (int i=1; i<=n; i++) {
				int[] s = sscoringMatrix4[s1e[i-1]];
				for (int j=1; j<=m; j++) {
					sI[i][j] = Math.max(sI[i-1][j]+gapExtend, sM[i-1][j]+gapFirst);
					sD[i][j] = Math.max(sD[i][j-1]+gapExtend, sM[i][j-1]+gapFirst);
					sM[i][j] = Math.max(sM[i-1][j-1]+s[s2e[j-1]], Math.max(sI[i][j], sD[i][j]));
				}
			}
			break;
		}
		return sM[n][m]/mult;
	}


	static void check(String[] args, int i)
	{
		if (args.length>i+1)
		{
			return;
		}
		System.err.println(args[i]+" takes an argument!");
		System.exit(-1);
	}

	public static void main(String[] args) throws IOException {

		int gap_open = -12;
		int gap_extend = -1;
		String scoring = "dayhoff.mat";
		HashMap<String,int[]> seqMapEnc = new HashMap<String, int[]>();
		HashMap<String,String> seqMap = new HashMap<String, String>();
		File pairs = new File("cathscop.inpairs");
		File seq = new File("domains.seqlib");

		boolean init = false;
		boolean intari = false;
		int mode = 4;
		
		for (int i=0; i<args.length; i++)
		{
			if (args[i].equals("-h") || args[i].equals("--help"))
			{
				usage();
				System.exit(0);
			}
			if (args[i].equals("-m"))
			{

				check(args,i);
				scoring = args[++i];
				continue;
			}
			if (args[i].equals("--go"))
			{
				check(args,i);
				gap_open = Integer.parseInt(args[++i]);
				continue;
			}
			if (args[i].equals("--reinit"))
			{
				check(args,i);
				init = args[++i].equals("1");
				continue;
			}
			if (args[i].equals("--mode"))
			{
				check(args, i);
				mode = Integer.parseInt(args[++i]);
				continue;
			}
			if (args[i].equals("--int"))
			{
				check(args, i);
				intari = args[++i].equals("1");
				continue;
			}
			if (args[i].equals("--ge"))
			{
				check(args,i);
				gap_extend = Integer.parseInt(args[++i]);
				continue;
			}

			if (args[i].equals("--pairs"))
			{
				check(args,i);
				pairs = new java.io.File(args[++i]);
				if (!pairs.exists())
				{
					System.err.println("pairfile: "+ pairs.getName()+" does not exist!");
					System.exit(-1);
				}
				continue;
			}
			if (args[i].equals("--seqlib"))
			{
				check(args,i);
				seq = new java.io.File(args[++i]);
			}
		}

		if (seq==null)
		{
			System.out.println("seqlib is not specified!");
			usage();
			System.exit(0);
		}
		if (pairs==null)
		{
			usage();
			System.exit(0);
		}


		ScoringMatrix m = ScoringMatrix.fromFile(new File(scoring));
		int maxLength = getSeqMap(seq,m,seqMap,seqMapEnc);
		Gotoh fg = new Gotoh(maxLength,m,gap_open,gap_extend);

		if (pairs!=null) {
			BufferedReader br = new BufferedReader(new FileReader(pairs));
			String line;
			String s1 = null; 
			String s2 = null; 
			int[] s1e = null; 
			int[] s2e = null; 
			double score;
			while ((line=br.readLine())!=null) {
				String id1 = line.substring(0,line.indexOf(' '));
				String id2 = line.substring(id1.length()+1);
				s1 = seqMap.get(id1);
				s2 = seqMap.get(id2);
				s1e = seqMapEnc.get(id1);
				s2e = seqMapEnc.get(id2);
				if (intari)
					score = fg.alignScaled(s1,s2,s1e, s2e, mode,init);
				else
					score = fg.align(s1,s2,s1e, s2e,mode, init);
				System.out.printf(Locale.US,"%s %s %.4f\n",id1,id2,score);
			}
		}

	}


	private static int[] getSequence(String ali, ScoringMatrix m) {
		int p = ali.indexOf(':');
		ali = ali.substring(p+1).trim().replace("-", "");
		return m.encodeSequence(ali);
	}

	private static String getId(String ali) {
		int p = ali.indexOf(':');
		if (p==-1) {
			System.err.println("alifile: not an alignment line "+ali);
			System.exit(-1);
		}
		return ali.substring(0,p);
	}

	static void usage() {
		System.out.println("usage for many pairs : --pairs <pairfile> --seqlib <seqlibfile> <optinal params>");
		System.out.println();
		System.out.println("	where seqlib contains lines with the format id:seq");
		System.out.println();
		System.out.println("		         queryid, targetid and the two first colums of pairfile has to be defined in the seqlibfile");
		System.out.println();
		System.out.println();
		System.out.println("        the optional parameters are:");
		System.out.println();
		System.out.println();
		System.out.println("                -m matrixfile in quasar format, default dayhoff.mat");
		System.out.println("                --go gapopen (default -12)");
		System.out.println("                --ge gapextend (default -1)");
		System.out.println();

	}

	private static int getSeqMap(File file, ScoringMatrix m, HashMap<String,String> re,HashMap<String,int[]> red) throws IOException {
		int max = 0;
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while ((line=br.readLine())!=null) {
			String[] p = line.split(":");
			if(p.length==2) {
				int[] s = m.encodeSequence(p[1]);
				red.put(p[0],s);
				re.put(p[0],p[1]);
				max = Math.max(max,s.length);
			}
		}
		br.close();
		return max;
	}

}


