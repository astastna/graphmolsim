package util;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Dictionary;
import java.util.HashMap;
import java.util.Locale;
import java.util.Properties;


public class ResultPrinter {

	/**
	 * where the results are printed
	 */
	private String resultFolder;


	/**
	 * the decimalformat for the editdistances found
	 */
	private DecimalFormat decFormat;
	
	/**
	 * the properties defined via GUI or properties file
	 */
	private Properties properties;

	
	/**
	 * constructs a ResultPrinter 
	 * @param resultFolder
	 * @param properties
	 */
	public ResultPrinter(String resultFolder, Properties properties) {
		this.resultFolder = resultFolder;
		this.properties = properties;
		
	}

	/**prints out the properties and the distance-matrix @param d 
	 * in the resultFolder as a json file
	 * 
	 * @author Aneta Šťastná
	 * @param matchingTime 
	 * @param numOfFails 
	 */
	public void printResult(double[][] d, long matchingTime, int numOfFails, String[] sourceNames, String[] targetNames, HashMap<Integer, String> indexSymbol, double[][] atomCostsMatrix) {
		this.decFormat = (DecimalFormat) NumberFormat
				.getInstance(Locale.ENGLISH);
		this.decFormat.applyPattern("0.00000");
		
		String resultName = this.resultFolder+ this.properties.getProperty("matching") + "-" + this.properties.getProperty("alpha") +".json";
		System.out.println(resultName);
		PrintWriter out;
		try {
			out = new PrintWriter(new FileOutputStream(resultName));
			out.println("{");
			out.println("\"properties\": {");
			out.println(" \"source\": {\n  \"file\": \"" + this.properties.getProperty("source")+"\",\n  \"graphs\": " +d.length+" \n },");
			out.println(" \"target\": {\n  \"file\": \"" + this.properties.getProperty("target")+"\",\n  \"graphs\": " +d[0].length+" \n },");
			//out.println("#Target graph set:\t"+this.properties.getProperty("target")+" ("+d[0].length+" graphs)");
			out.println(" \"algorithm\": \""+this.properties.getProperty("matching")+ "\",");
			//out.println("#Graph edit distance procedure:\t"+this.properties.getProperty("matching"));
			if (this.properties.getProperty("matching").equals("Beam")){
				out.println(" \"s\": "+ Integer.parseInt(properties.getProperty("s")) + ",");
			}
			int undirected = Integer
					.parseInt(properties.getProperty("undirected"));
			out.print(" \"edges\": ");
			if (undirected == 1){
				out.println("\"undirected\",");
			} else {
				out.println("\"directed\",");;
			}
			out.println(" \"node delete/insert\": "+this.properties.getProperty("node") + ",");
			out.println(" \"edge delete/insert\": "+this.properties.getProperty("edge")+",");
			out.println(" \"node-edge weight ratio\": "+this.properties.getProperty("alpha") + ",");
			int numOfNodeAttr = Integer.parseInt(properties
					.getProperty("numOfNodeAttr"));
			int numOfEdgeAttr = Integer.parseInt(properties
					.getProperty("numOfEdgeAttr"));
			
			if (numOfNodeAttr > 0) out.println(" \"node attributes\": [");
			for (int i = 0; i < numOfNodeAttr; i++) {
				out.println("  {\n  \"name\": \""+ properties.getProperty("nodeAttr" + i)+"\",");
				out.println("  \"cost function\": \""+properties.getProperty("nodeCostType" + i)+"\",");
				if (properties.getProperty("nodeCostType" + i).equals("discrete")){
					out.println("  \"mu\":  \""+properties.getProperty("nodeCostMu" + i)+"\"\n  \"nu\": \""+properties.getProperty("nodeCostNu" + i)+"\",");
				}
				out.println("  \"attr importance\": "+properties.getProperty("nodeAttr" + i + "Importance") + "");
				if (i+1 == numOfNodeAttr) out.println("  }");
				else out.println("  },");
			}
			if (numOfNodeAttr > 0) out.println(" ],");
			/*if (numOfNodeAttr==0){
				out.println("#No attributes for nodes defined");
			}*/

			if (numOfEdgeAttr > 0) out.println(" \"edge attributes\":[");
			for (int i = 0; i < numOfEdgeAttr; i++) {
				out.println("  {");
				out.println("   \"name\": \""+ properties.getProperty("edgeAttr" + i)+"\",");
				out.println("   \"cost Function\": \""+properties.getProperty("edgeCostType" + i)+"\",");
				out.println("   \"attr importance\": "+properties.getProperty("edgeAttr" + i + "Importance"));
				if (i+1 == numOfEdgeAttr) out.println("  }");
				else out.println("  },");
			}
			/*if (numOfEdgeAttr==0){
				out.println("#No attributes for edges defined");
			}*/
			if (numOfEdgeAttr > 0) out.println("  ],");

			double squareRootNodeCosts = Double.parseDouble(properties
					.getProperty("pNode"));
			int multiplyNodeCosts = Integer.parseInt(properties
					.getProperty("multiplyNodeCosts"));
			double squareRootEdgeCosts = Double.parseDouble(properties
					.getProperty("pEdge"));
			int multiplyEdgeCosts = Integer.parseInt(properties
					.getProperty("multiplyEdgeCosts"));
			
			if (multiplyNodeCosts==1){
				out.println(" \"node costs count\" : \"multiply\",");
			} else {
				out.println(" \"node costs count\" : \"add\",");
			}
			if (multiplyEdgeCosts==1){
				out.println(" \"edge costs count\" :  \"multiply\",");
			} else {
				out.println(" \"edge costs count\" : \"add\",");
			}	
			out.println(" \"square root node costs\" : "+squareRootNodeCosts+",");
		
			out.println(" \"square root edge costs\" : "+squareRootEdgeCosts+",");
			
			int simKernel=Integer.parseInt(properties.getProperty("simKernel"));
			
			out.println(" \"time\" : "+ matchingTime + ",");
			out.println(" \"fails\" : "+ numOfFails + ",");
			
			out.println("  \"element change cost matrix\": {");
			boolean firstLine = true;
			for (int i = 0; i < atomCostsMatrix.length; i++){
				
				if (firstLine){ //dealing with the first line
					out.print("   \"elements\": [ ");
					for (int j = 0; j < atomCostsMatrix[0].length; j++){
						if ((j+1) == atomCostsMatrix[0].length){
							out.print( "\"" + indexSymbol.get(j) + "\"");//last element
						} else {
							out.print( "\"" + indexSymbol.get(j) + "\", ");
						}
					}
					out.println("],");
					firstLine = false;
				}
				
				// printing lines with source headers
				out.print("   \""+ indexSymbol.get(i) +"\": [ "); 
				for (int j = 0; j < atomCostsMatrix[0].length; j++){
					if (j+1 == atomCostsMatrix[0].length){//last element
						out.print(decFormat.format(atomCostsMatrix[i][j]));
					} else{
					out.print(decFormat.format(atomCostsMatrix[i][j])+", ");
					}
				}
				if ( (i+1) == atomCostsMatrix.length){
					out.println("]");
				} else {
					out.println("], ");
				}
			}
			out.println("  }");
			
			out.println("},");//end of properties
			
			out.println(" \"results\": {");
			out.print("  \"type\" : ");
			switch (simKernel){
				case 0:
					out.println("\"distance matrix\",");
					break;
				case 1:
					out.println("\"similarity matrix (-d^2)\",");break;
				case 2:
					out.println("\"similarity matrix (-d)\",");break;
				case 3:
					out.println("\"similarity matrix tanh(-d)\",");break;
				case 4:
					out.println("\"similarity matrix exp(-d)\",");break;
			}
			
			
			/* distance matrix in json format:
			 * {
			 * "targets": [target1, target2, ... ],
			 * "source1": [val1, val2, ... ],
			 * .
			 * .
			 * .
			 * }
			 * 
			 * @author Aneta Šťastná
			 * 
			 */
			out.println("  \"matrix\": {");
			boolean firstLine1 = true;
			for (int i = 0; i < d.length; i++){
				
				if (firstLine1){ //dealing with the first line
					out.print("   \"targets\": [ ");
					for (int j = 0; j < d[0].length; j++){
						if ((j+1) == d[0].length){
							out.print( "\"" + targetNames[j] + "\"");//last element
						} else {
							out.print( "\"" + targetNames[j] + "\", ");
						}
					}
					out.println("],");
					firstLine1 = false;
				}
				
				// printing lines with source headers
				out.print("   \""+ sourceNames[i] +"\": [ "); 
				for (int j = 0; j < d[0].length; j++){
					if (j+1 == d[0].length){//last element
						out.print(decFormat.format(d[i][j]));
					} else{
					out.print(decFormat.format(d[i][j])+", ");
					}
				}
				if ( (i+1) == d.length){
					out.println("]");
				} else {
					out.println("], ");
				}
			}
			out.println("  }");
			
			out.println(" }"); //end of results
			out.println("}"); //end of document
			
			out.flush();
			out.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}

	public void printResult(double[][] d,
			double[][] a, long matchingTime, int numOfFails) {
		this.decFormat = (DecimalFormat) NumberFormat
				.getInstance(Locale.ENGLISH);
		this.decFormat.applyPattern("0.00000");
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
		Calendar cal = Calendar.getInstance(); // +"_"+dateFormat.format(cal.getTime())
		String resultName = this.resultFolder+".txt";
		System.out.println(resultName);
		PrintWriter out;
		try {
			out = new PrintWriter(new FileOutputStream(resultName));
			out.println("#*** The properties of the matching ***");
			out.println("#Source graph set:\t"+this.properties.getProperty("source")+" ("+d.length+" graphs)");
			out.println("#Target graph set:\t"+this.properties.getProperty("target")+" ("+d[0].length+" graphs)");
			out.println("#Graph edit distance procedure:\t"+this.properties.getProperty("matching"));
			if (this.properties.getProperty("matching").equals("Beam")){
				out.println("#s = "+ Integer.parseInt(properties.getProperty("s")));
			}
			int undirected = Integer
					.parseInt(properties.getProperty("undirected"));
			out.print("#Edge mode:\t\t");
			if (undirected == 1){
				out.println("undirected");
			} else {
				out.println("directed");;
			}
			out.println("#Cost for node deletion/insertion:\t"+this.properties.getProperty("node"));
			out.println("#Cost for edge deletion/insertion:\t"+this.properties.getProperty("edge")+"");
			out.println("#Alpha weighting factor between node and edge costs:\t"+this.properties.getProperty("alpha"));
			int numOfNodeAttr = Integer.parseInt(properties
					.getProperty("numOfNodeAttr"));
			int numOfEdgeAttr = Integer.parseInt(properties
					.getProperty("numOfEdgeAttr"));
			
			
			for (int i = 0; i < numOfNodeAttr; i++) {
				out.print("#Node attribute "+i+":\t"+ properties.getProperty("nodeAttr" + i)+";\t");
				out.print("Cost function:\t"+properties.getProperty("nodeCostType" + i)+";\t");
				if (properties.getProperty("nodeCostType" + i).equals("discrete")){
					out.print("mu = "+properties.getProperty("nodeCostMu" + i)+" nu = "+properties.getProperty("nodeCostNu" + i)+";\t");
				}
				out.println("Soft factor:\t"+properties.getProperty("nodeAttr" + i + "Importance"));
			}
			if (numOfNodeAttr==0){
				out.println("#No attributes for nodes defined");
			}

			for (int i = 0; i < numOfEdgeAttr; i++) {
				out.print("#Edge Attribute "+i+":\t"+ properties.getProperty("edgeAttr" + i)+";\t");
				out.print("Cost Function:\t"+properties.getProperty("edgeCostType" + i)+";\t");
				out.println("Soft Factor:\t"+properties.getProperty("edgeAttr" + i + "Importance"));
			}
			if (numOfEdgeAttr==0){
				out.println("#No attributes for edges defined");
			}

			double squareRootNodeCosts = Double.parseDouble(properties
					.getProperty("pNode"));
			int multiplyNodeCosts = Integer.parseInt(properties
					.getProperty("multiplyNodeCosts"));
			double squareRootEdgeCosts = Double.parseDouble(properties
					.getProperty("pEdge"));
			int multiplyEdgeCosts = Integer.parseInt(properties
					.getProperty("multiplyEdgeCosts"));
			
			if (multiplyNodeCosts==1){
				out.println("#Individual node costs are multiplied");
			} else {
				out.println("#Individual node costs are added");
			}
			if (multiplyEdgeCosts==1){
				out.println("#Individual edge costs are multiplied");
			} else {
				out.println("#Individual edge costs are added");
			}	
			out.println("#(Combined node cost)^(1/"+squareRootNodeCosts+")");
		
			out.println("#(Combined edge cost)^(1/"+squareRootEdgeCosts+")");
			
			int simKernel=Integer.parseInt(properties.getProperty("simKernel"));
			
			out.println("#Complete Matching Time: "+ matchingTime);
			out.println("#Number of Fails: "+ numOfFails);
			
			
			
			for (int i = 0; i < d.length; i++){
				for (int j = 0; j < d[0].length; j++){
//					out.print(decFormat.format(d[i][j])+" ");
//					out.print(decFormat.format(a[i][j])+"; ");
//					if (j+1 < d[0].length){
//						out.print(",");
//					}
				}
				out.println();
			}
			
			out.flush();
			out.close();
			//this function is not being used!!!
			
			
		} catch (FileNotFoundException e) {
			System.out.println("Output file not found.");
			e.printStackTrace();
		}
		
	}

}
