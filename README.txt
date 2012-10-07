A demo of elastic image warping in HTML5 using JavaScript and the Canvas API. 2D
image warping code using deformation meshes adopted to JS from Java code pub-
lished by Jerry Huxtable at http://www.jhlabs.com/ip/index.html. Mesh construct-
ion and deformation code and 1D image stretching code (embedded in the HTML file
rather than in a JS file) written from scratch. 

Demo image taken from Wikipedia / Wikimedia Commons: 
http://en.wikipedia.org/wiki/File:African_landscape.jpg

Due to the security model, this demo will not work in most browsers if you 
simply open the HTML from disk. Both the warp.html page and the example JPEG 
image MUST be hosted on the same server.
For a easy solution, you could run e.g. "python -m SimpleHTTPServer" in the 
./src directory and access the page at http://localhost:8000/warp.html
