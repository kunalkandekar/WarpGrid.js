//image warping stuff
/*
    Adopted from Java code by Jerry Huxtable at http://www.jhlabs.com/ip/index.html
*/
var m00 = -0.5;
var m01 =  1.5;
var m02 = -1.5;
var m03 =  0.5;
var m10 =  1.0;
var m11 = -2.5;
var m12 =  2.0;
var m13 = -0.5;
var m20 = -0.5;
var m22 =  0.5;
var m31 =  1.0;
	
function setval(arr, val) {
	for(var i = 0, n = arr.length; i < n; i++) { arr[i] = val; }
}

//use float32/int32 arrays if browser has them
var FloatArray = Array;
if(typeof Float32Array == 'function') { 
    FloatArray = Float32Array;
}
var IntArray = Array;
if(typeof Int32Array == 'function') { 
    IntArray = Int32Array;
}

function newFloatArray(size, defval) {
   var arr = new FloatArray(size);//Array(size);
   if(defval) {
       setval(arr, defval);
   }
   return arr;
}

function newIntArray(size, defval) {
   var arr = new Int32Array(size);//Array(size);
   if(defval) {
       setval(arr, defval);
   }
   return arr;
}	
    
function packToIntArray(unpacked, packed) {
    if(!packed) {
        packed = newIntArray(unpacked.length / 4);
    }
    for(var i = 0, n = packed.length; i < n; i++) {
        var j = i * 4;
        packed[i] = (unpacked[j + 3] << 24) |
					(unpacked[j] << 16) |   
					(unpacked[j + 1] << 8) |
					(unpacked[j + 2]);
    } 
    //debug("Packed ["+unpacked[0]+","+unpacked[1]+","+unpacked[2]+","+unpacked[2]+"] to "+packed[0]); 
    return packed;
}

function unpackFromIntArray(packed, unpacked) {
    if(!unpacked) {
        unpacked = newIntArray(packed.length * 4);
    }
    var index = 0;
    for(var i = 0, n = packed.length; i < n; i++) {
        var rgb = packed[i];
		unpacked[index++] = (rgb >> 16) & 0xff;
		unpacked[index++] = (rgb >> 8) & 0xff;
		unpacked[index++] = (rgb & 0xff);
		unpacked[index++] = (rgb >> 24) & 0xff;
    }
    //debug("Unpacked: " + packed[0]+" to  ["+unpacked[0]+","+unpacked[1]+","+unpacked[2]+","+unpacked[2]+"]"); 
    return unpacked;
}

//TODO: Figure out how to port this to WebGL
function resample(source, dest, length, offset, stride, out) {
	//public static void resample(int[] source, int[] dest, int length, int offset, int stride, float[] out)
	var i = 0, j = 0;
	var sizfac = 0.0;
	var inSegment = 0.0;
	var outSegment = 0.0;
	var a = 0, r = 0, g = 0, b = 0, nextA = 0, nextR = 0, nextG = 0, nextB = 0;
	var aSum = 0.0, rSum = 0.0, gSum = 0.0, bSum = 0.0;
	//float[] in;
	var srcIndex = offset;
	var destIndex = offset;
	var lastIndex = source.length;
	var rgb = 0;

	var inp = newFloatArray(length+2);
	i = 0;
	for (j = 0; j < length; j++) {
		while (out[i+1] < j)
			i++;
		inp[j] = i + 1.0 * (j - out[i]) / (out[i + 1] - out[i]);
        //inp[j] = clamp( inp[j], 0, length-1 );
	}
	inp[length] = length;
	inp[length+1] = length;

	inSegment  = 1.0;
	outSegment = inp[1];
	sizfac = outSegment;
	aSum = rSum = gSum = bSum = 0.0;
	rgb = source[srcIndex];
	a = (rgb >> 24) & 0xff;
	r = (rgb >> 16) & 0xff;
	g = (rgb >> 8) & 0xff;
	b = rgb & 0xff;
	srcIndex += stride;
	rgb = source[srcIndex];
	nextA = (rgb >> 24) & 0xff;
	nextR = (rgb >> 16) & 0xff;
	nextG = (rgb >> 8) & 0xff;
	nextB = rgb & 0xff;
	srcIndex += stride;
	i = 1;

	while (i <= length) {
		var aIntensity = inSegment * a + (1.0 - inSegment) * nextA;
		var rIntensity = inSegment * r + (1.0 - inSegment) * nextR;
		var gIntensity = inSegment * g + (1.0 - inSegment) * nextG;
		var bIntensity = inSegment * b + (1.0 - inSegment) * nextB;
		if (inSegment < outSegment) {
			aSum += (aIntensity * inSegment);
			rSum += (rIntensity * inSegment);
			gSum += (gIntensity * inSegment);
			bSum += (bIntensity * inSegment);
			outSegment -= inSegment;
			inSegment = 1.0;
			a = nextA;
			r = nextR;
			g = nextG;
			b = nextB;
			if (srcIndex < lastIndex)
				rgb = source[srcIndex];
			nextA = (rgb >> 24) & 0xff;
			nextR = (rgb >> 16) & 0xff;
			nextG = (rgb >> 8) & 0xff;
			nextB = rgb & 0xff;
			srcIndex += stride;
		} else {
			aSum += (aIntensity * outSegment);
			rSum += (rIntensity * outSegment);
			gSum += (gIntensity * outSegment);
			bSum += (bIntensity * outSegment);
			dest[destIndex] = 
				(Math.round(Math.min(aSum/sizfac, 255) << 24)) |
				(Math.round(Math.min(rSum/sizfac, 255) << 16)) |
				(Math.round(Math.min(gSum/sizfac, 255) << 8)) |
				Math.round(Math.min(bSum/sizfac, 255));
			destIndex += stride;
			aSum = rSum = gSum = bSum = 0.0;
			inSegment -= outSegment;
			outSegment = inp[i+1] - inp[i];
			sizfac = outSegment;
			i++;
		}
	}
}

function WarpGrid(rows, cols, w, h) {        	        
   //public WarpGrid(int rows, int cols, int w, int h) {
	this.rows = rows;
	this.cols = cols;
	this.xGrid = newFloatArray(rows*cols);
	this.yGrid = newFloatArray(rows*cols);
	var index = 0;
	for (var row = 0; row < rows; row++) {
		for (var col = 0; col < cols; col++) {
			this.xGrid[index] = 1.0 * (col*(w-1)/(cols-1));
			this.yGrid[index] = 1.0 * (row*(h-1)/(rows-1));
			index++;
		}
	}
}

WarpGrid.prototype.warp = function(inPixels, cols, rows, sourceGrid, destGrid, outPixels) {
   //(int[] inPixels, int cols, int rows, WarpGrid sourceGrid, WarpGrid destGrid, int[] outPixels) {
    /*try { */
		var x = 0, y = 0; //int
		var u = 0, v = 0; //int
		//var intermediate = new Array();
		//WarpGrid splines;

		if (sourceGrid.rows != destGrid.rows || sourceGrid.cols != destGrid.cols)
			throw ("source and destination grids are different sizes");

		var size = Math.max(cols, rows);
		var xrow = newFloatArray(size);
		var yrow = newFloatArray(size);
		var scale  = newFloatArray(size + 1);
		var interpolated = newFloatArray(size + 1);
		//init to 0
		setval(interpolated, 0);

		var gridCols = sourceGrid.cols;
		var gridRows = sourceGrid.rows;

		splines = new WarpGrid(rows, gridCols, 1, 1);

		for (u = 0; u < gridCols;u++) {
			var i = u;

			for (v = 0; v < gridRows;v++) {
				xrow[v] = sourceGrid.xGrid[i];
				yrow[v] = sourceGrid.yGrid[i];
				i += gridCols;
			}

			this.interpolateSpline(yrow, xrow, 0, gridRows, interpolated, 0, rows);

			i = u;
			for (y = 0;y < rows;y++) {
				splines.xGrid[i] = interpolated[y];
				i += gridCols;
			}
		}

		for (u = 0; u < gridCols;u++) {
			var i = u;

			for (v = 0; v < gridRows;v++) {
				xrow[v] = destGrid.xGrid[i];
				yrow[v] = destGrid.yGrid[i];
				i += gridCols;
			}

			this.interpolateSpline(yrow, xrow, 0, gridRows, interpolated, 0, rows);

			i = u;
			for (y = 0;y < rows; y++) {
				splines.yGrid[i] = interpolated[y];
				i += gridCols;
			}
		}

		/* first pass: warp x using splines */
		intermediate = newIntArray(rows*cols);  //int

		var offset = 0;
	    timeIt("warp: pre-resample 1 = "+rows);
		for (y = 0; y < rows; y++) {
			/* fit spline to x-intercepts;resample over all cols */
			this.interpolateSpline(splines.xGrid, splines.yGrid, offset, gridCols, scale, 0, cols);
			scale[cols] = cols;
			resample(inPixels, intermediate, cols, y*cols, 1, scale);
			offset += gridCols;
		}
	    timeIt("warp: resample 1");

		/* create table of y-intercepts for intermediate mesh's hor splines */

		splines = new WarpGrid(gridRows, cols, 1, 1);

		offset = 0;
		var offset2 = 0;
		for (v = 0; v < gridRows; v++) {
			this.interpolateSpline(sourceGrid.xGrid, sourceGrid.yGrid, offset, gridCols, splines.xGrid, offset2, cols);
			offset += gridCols;
			offset2 += cols;
		}

		offset = 0;
		offset2 = 0;
		for (v = 0; v < gridRows; v++) {
			this.interpolateSpline(destGrid.xGrid, destGrid.yGrid, offset, gridCols, splines.yGrid, offset2, cols);
			offset += gridCols;
			offset2 += cols;
		}
        	        
		/* second pass: warp y */

		for (x = 0; x < cols; x++) {
			var i = x;
			
			for (v = 0; v < gridRows; v++) {
				xrow[v] = splines.xGrid[i];
				yrow[v] = splines.yGrid[i];
				i += cols;
			}

			this.interpolateSpline(xrow, yrow, 0, gridRows, scale, 0, rows);
			scale[rows] = rows;
			resample(intermediate, outPixels, rows, x, cols, scale);
		}
        	timeIt("warp: 2ndpass");
    /*}
    catch (e) {
    	debug(e);
    }*/
}

WarpGrid.prototype.interpolateSpline = function(xKnots, yKnots, offset, length, splineY, splineOffset, splineLength) {
//protected void interpolateSpline(float[] xKnots, float[] yKnots, int offset, int length, float[] splineY, int splineOffset, int splineLength) {
	var index = offset;
	var end = offset+length-1;
	var x0 = 0.0, x1 = 0.0;
	var k0 = 0.0, k1 = 0.0, k2 = 0.0, k3 = 0.0;
	var c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0;

	x0 = xKnots[index];
	k0 = k1 = k2 = yKnots[index];
	x1 = xKnots[index+1];
	k3 = yKnots[index+1];

	for (var i = 0;i < splineLength;i++) {
		if (index <= end && i > xKnots[index]) {
			k0 = k1;
			k1 = k2;
			k2 = k3;
			x0 = xKnots[index];
			index++;
			if ( index <= end )
				x1 = xKnots[index];
			if ( index < end )
				k3 = yKnots[index+1];
			else
				k3 = k2;
		}
		t = 1.0 * ((i - x0) / (x1 - x0));
		c3 = m00*k0 + m01*k1 + m02*k2 + m03*k3;
		c2 = m10*k0 + m11*k1 + m12*k2 + m13*k3;
		c1 = m20*k0 + m22*k2;
		c0 = m31*k1;
		
		splineY[splineOffset+i] = ((c3*t + c2)*t + c1)*t + c0;
	}
}