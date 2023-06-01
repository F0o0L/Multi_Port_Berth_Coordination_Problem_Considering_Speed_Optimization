package CLPSO_MVNS;


public class run {
	public static void main(String[] args) {
		long start = System.currentTimeMillis(); 
		String fileName="F:\\JavaWorkplace\\MPBAPJavaVersion\\src\\main\\java\\CLPSO_MVNS\\ini.txt";
		CLPSO_MVNS clpso_mvns=new CLPSO_MVNS(500,500,fileName);
//		CLPSO_MVNS clpso_mvns=new CLPSO_MVNS(500,500,3,3,5);
		clpso_mvns.run();
		long end = System.currentTimeMillis(); 
		System.out.println(end-start);
	}
}
