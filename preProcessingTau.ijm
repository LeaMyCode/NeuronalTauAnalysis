close("*");

//Open folders (first for the images you want to analyse, second to save the treshold images and the excel files with the results; you have to save this in a new ordner) 
dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Select destination for saving the analysis data");
list = getFileList(dir1);
print(list.length);

//set Measurement settings
run("Set Measurements...", "area centroid perimeter redirect=None decimal=3");

for (j=0; j<list.length; j++) {
	showProgress(j+1, list.length);
	if(list[j].matches('.*.zip') || list[j].matches('.*.tif') || list[j].matches('.*.czi')) {
		picture = list[j];
		open(""+dir1+picture+"");
	
	//get title for the image
		title = getTitle();
	run("Z Project...", "projection=[Max Intensity]");
	Stack.getDimensions(width, height, channels, slices, frames); 
	if (channels == 2) {
	//split channels 
		run("Split Channels");
		selectWindow ("C1-MAX_" + picture);
		rename("AT8");
		selectWindow ("C2-MAX_" + picture);
		rename("Neuron");
		
		//Setting measurments
		run("Set Measurements...", "area mean min perimeter integrated redirect=None decimal=3");
		
		//Tresholding of the neuron
		selectWindow("Neuron"); 
		run("Duplicate...", "title=Neuron-1");
		selectWindow("Neuron-1"); 
		run("Enhance Contrast", "saturated=0.35");
		run("Despeckle");
		setAutoThreshold("Li dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		
		waitForUser("Please delete all background from your selected neuron.");
		
		//run("Fill Holes");
		
		selectWindow("Neuron-1"); 
		run("Duplicate...", "title=Neuron_threshold");
		selectWindow("Neuron_threshold");
		saveAs("PNG", dir2 + title + "-Neuron_threshold.png");
		run("Close");
		
		//Getting information about the soma
		selectWindow("Neuron-1");
		setTool("freehand");
		waitForUser("Please mark the SOMA of your cell, add the ROI to your Roi manager and then press ok.");
	    for (r = 0; r < roiManager("count"); r++) {
	        roiManager("select", r);
	        run("Measure");
	    }
	    roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-soma_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-soma_ROI.csv");
		run("Close");
		
		//cut the soma out for the dendrite analysis
		roiManager("Select", 0);
		setBackgroundColor(255, 255, 255);
		run("Clear", "slice");
		
		
		roiManager("reset");
		run("Select None");
		
		//Getting information about the dendrites
		selectWindow("Neuron-1");
		waitForUser("Make sure the soma is cut out.");
		selectWindow("Neuron-1");
		run("Analyze Particles...", "size=30-Infinity add");
		print(roiManager("count"));

	    for (r = 0; r < roiManager("count"); r++) {
	        roiManager("select", r);
	        run("Measure");
	    }
	    roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-dendrites_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-dendrites_ROI.csv");
		run("Close");
		roiManager("reset");
		
		
		//Getting information about the nucelus
		selectWindow("Neuron");
		setTool("freehand");
		waitForUser("Please mark the NUCLEUS of your cell, add the ROI to your Roi manager and then press ok.");
		roiManager('select', "*");
		run("Measure");
		roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-nucleus_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-nucleus_ROI.csv");
		run("Close");
		roiManager("reset");
		

	//hyperphosTau processing
	//copy
		selectWindow("AT8");
		run("Enhance Contrast", "saturated=0.75");
		run("Duplicate...", "title=AT8-1");
	
	//unsharp mask
		selectWindow("AT8");
		run("Unsharp Mask...", "radius=3 mask=0.60");
		
	//Gaussian Blur
		selectWindow("AT8-1");
		run("Gaussian Blur...", "sigma=25");
	//Image calculator
		imageCalculator("Subtract create stack","AT8", "AT8-1");
			
	//Thresholding		
		selectWindow("Result of AT8");
		setAutoThreshold("Moments dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Despeckle");
		run("Watershed");
	//Analyze particles
		run("Analyze Particles...", "size= 0-infinity show=Nothing summarize add");
		
	//Measure particle size and count particles
		selectWindow("Summary");
		run("Close");
		
		for (r = 0; r < roiManager("count"); r++) {
		        roiManager("select", r);
		        run("Measure");
		    }
		
	//Save
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-AT8.csv");
		run("Close");
		selectWindow("Result of AT8");
		saveAs("PNG", dir2 + title + "-AT8_treshold.png");
		roiManager('select', "*");
		roiManager("Save", dir2 + title+ "AT8_ROI.zip");
		
	//cleanup roimanager
	roiManager("reset");


		}
	
	if (channels == 3) {
	//split channels 
		run("Split Channels");
		selectWindow ("C1-MAX_" + picture);
		rename("totalTau");
		selectWindow ("C2-MAX_" + picture);
		rename("AT8");
		selectWindow ("C3-MAX_" + picture);
		rename("Neuron");
		
		//Setting measurments
		run("Set Measurements...", "area mean min perimeter integrated redirect=None decimal=3");
		
		//Tresholding of the neuron
		selectWindow("Neuron"); 
		run("Duplicate...", "title=Neuron-1");
		selectWindow("Neuron-1"); 
		run("Enhance Contrast", "saturated=0.35");
		run("Despeckle");
		setAutoThreshold("Li dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		
		waitForUser("Please delete all background from your selected neuron.");
		
		//run("Fill Holes");
		
		selectWindow("Neuron-1"); 
		run("Duplicate...", "title=Neuron_threshold");
		selectWindow("Neuron_threshold");
		saveAs("PNG", dir2 + title + "-Neuron_threshold.png");
		run("Close");
		
		//Getting information about the soma
		selectWindow("Neuron-1");
		setTool("freehand");
		waitForUser("Please mark the SOMA of your cell, add the ROI to your Roi manager and then press ok.");
	    for (r = 0; r < roiManager("count"); r++) {
	        roiManager("select", r);
	        run("Measure");
	    }
	    roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-soma_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-soma_ROI.csv");
		run("Close");
		
		//cut the soma out for the dendrite analysis
		roiManager("Select", 0);
		setBackgroundColor(255, 255, 255);
		run("Clear", "slice");
		
		
		roiManager("reset");
		run("Select None");
		
		//Getting information about the dendrites
		selectWindow("Neuron-1");
		waitForUser("Make sure the soma is cut out.");
		selectWindow("Neuron-1");
		run("Analyze Particles...", "size=30-Infinity add");
		print(roiManager("count"));

	    for (r = 0; r < roiManager("count"); r++) {
	        roiManager("select", r);
	        run("Measure");
	    }
	    roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-dendrites_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-dendrites_ROI.csv");
		run("Close");
		roiManager("reset");
		
		
		//Getting information about the nucelus
		selectWindow("Neuron");
		setTool("freehand");
		waitForUser("Please mark the NUCLEUS of your cell, add the ROI to your Roi manager and then press ok.");
		roiManager('select', "*");
		run("Measure");
		roiManager('select', "*");
		roiManager("Save", dir2 + title+ "-nucleus_ROI.zip");
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-nucleus_ROI.csv");
		run("Close");
		roiManager("reset");
		
	//total Tau processing
	//copy
		selectWindow("totalTau");
		run("Duplicate...", "title=totalTau-1");
	
	//unsharp mask
		selectWindow("totalTau");
		run("Unsharp Mask...", "radius=3 mask=0.60");
		
	//Gaussian Blur
		selectWindow("totalTau-1");
		run("Gaussian Blur...", "sigma=25");
	//Image calculator
		imageCalculator("Subtract create stack","totalTau", "totalTau-1");
			
	//Thresholding		
		selectWindow("Result of totalTau");
		setAutoThreshold("Moments dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Despeckle");
		run("Watershed");
	//Analyze particles
		run("Analyze Particles...", "size= 0-infinity show=Nothing summarize add");
		
	//Measure particle size and count particles
		selectWindow("Summary");
		run("Close");
		
		for (r = 0; r < roiManager("count"); r++) {
		        roiManager("select", r);
		        run("Measure");
		    }
		    
	//Save
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-totalTau.csv");
		run("Close");
		selectWindow("Result of totalTau");
		saveAs("PNG", dir2 + title + "-totalTau_treshold.png");
		roiManager('select', "*");
		roiManager("Save", dir2 + title+ "totalTau_ROI.zip");
		
	//cleanup results and roimanager
		roiManager("reset");
		run("Close");
	
	//hyperphosTau processing
	//copy
		selectWindow("AT8");
		run("Enhance Contrast", "saturated=0.75");
		run("Duplicate...", "title=AT8-1");
	
	//unsharp mask
		selectWindow("AT8");
		run("Unsharp Mask...", "radius=3 mask=0.60");
		
	//Gaussian Blur
		selectWindow("AT8-1");
		run("Gaussian Blur...", "sigma=25");
	//Image calculator
		imageCalculator("Subtract create stack","AT8", "AT8-1");
			
	//Thresholding		
		selectWindow("Result of AT8");
		setAutoThreshold("Moments dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Despeckle");
		run("Watershed");
	//Analyze particles
		run("Analyze Particles...", "size= 0-infinity show=Nothing summarize add");
		
	//Measure particle size and count particles
		selectWindow("Summary");
		run("Close");
		
		for (r = 0; r < roiManager("count"); r++) {
		        roiManager("select", r);
		        run("Measure");
		    }
		
	//Save
		selectWindow("Results");
		saveAs("Measurements", dir2 + title + "-AT8.csv");
		run("Close");
		selectWindow("Result of AT8");
		saveAs("PNG", dir2 + title + "-AT8_treshold.png");
		roiManager('select', "*");
		roiManager("Save", dir2 + title+ "AT8_ROI.zip");
		
	//cleanup roimanager
	roiManager("reset");


		}
	}

//Clean-up to prepare for next image
	run("Close All");
	close("*");

	if (isOpen("Log")) {
	     selectWindow("Log");
	     run("Close");
	}
	if (isOpen("Summary")) {
	     selectWindow("Summary");
	     run("Close");
	}
	if (isOpen("Results")) {
	     selectWindow("Results");
	     run("Results");
	}
	
  }
print("Jeah, finished!");

		


