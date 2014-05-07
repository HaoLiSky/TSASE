//V3.01.A - http://www.openjs.com/scripts/jx/
jx = {
	//Create a xmlHttpRequest object - this is the constructor. 
	getHTTPObject : function() {
		var http = false;
		//Use IE's ActiveX items to load the file.
		if(typeof ActiveXObject != 'undefined') {
			try {http = new ActiveXObject("Msxml2.XMLHTTP");}
			catch (e) {
				try {http = new ActiveXObject("Microsoft.XMLHTTP");}
				catch (E) {http = false;}
			}
		//If ActiveX is not available, use the XMLHttpRequest of Firefox/Mozilla etc. to load the document.
		} else if (window.XMLHttpRequest) {
			try {http = new XMLHttpRequest();}
			catch (e) {http = false;}
		}
		return http;
	},
	// This function is called from the user's script. 
	//Arguments - 
	//	url	- The url of the serverside script that is to be called. Append all the arguments to 
	//			this url - eg. 'get_data.php?id=5&car=benz'
	//	callback - Function that must be called once the data is ready.
	//	format - The return type for this function. Could be 'xml','json' or 'text'. If it is json, 
	//			the string will be 'eval'ed before returning it. Default:'text'
	load : function (url,callback,format) {
		var http = this.init(); //The XMLHttpRequest object is recreated at every call - to defeat Cache problem in IE
		if(!http||!url) return;
		if (http.overrideMimeType) http.overrideMimeType('text/xml');

		if(!format) var format = "text";//Default return type is 'text'
		format = format.toLowerCase();
		
		//Kill the Cache problem in IE.
		var now = "uid=" + new Date().getTime();
		url += (url.indexOf("?")+1) ? "&" : "?";
		url += now;

		http.open("GET", url, true);

		http.onreadystatechange = function () {//Call a function when the state changes.
			if (http.readyState == 4) {//Ready State will be 4 when the document is loaded.
				if(http.status == 200) {
					var result = "";
					if(http.responseText) result = http.responseText;
					
					//If the return is in JSON format, eval the result before returning it.
					if(format.charAt(0) == "j") {
						//\n's in JSON string, when evaluated will create errors in IE
						result = result.replace(/[\n\r]/g,"");
						result = eval('('+result+')'); 
					}
	
					//Give the data to the callback function.
					if(callback) callback(result);
				} else { //An error occured
					if(error) error(http.status);
				}
			}
		}
		http.send(null);
	},
	init : function() {return this.getHTTPObject();}
}


function xyz(canvas) {

    var self = this;
    
    self.atomColor = [0.8, 0.8, 1.0];
    self.clearColor = [1, 1, 1];
    self.animator = null;
    self.radius = 1.0;
    self.index = 0;
    self.fps = 0;
    self.mfps = 0;

    self.event_mousemove = function(e) {
        if(self.lastmousex == undefined) {
            self.lastmousex = e.layerX;
            self.lastmousey = e.layerY;
            return;
        }
        var dx = e.layerX - self.lastmousex;
        var dy = e.layerY - self.lastmousey;
        self.lastmousex = e.layerX;
        self.lastmousey = e.layerY;
        if(self.button1) {
            self.rotx(-dy * 0.02);
            self.roty(dx * 0.02);
            self.draw()
        }
    }

    self.event_mousedown = function(e) {
        self.canvas.focus();
        e.preventDefault();
        if(e.which == 1) {
            self.button1 = true;
        }
    }
    
    self.event_mouseup = function(e) {
        if(e.which == 1) {
            self.button1 = false;
        }
    }

    self.event_mousewheel = function(e) {
        e.preventDefault();
        if(e.wheelDelta > 0) {
            self.scale *= 1.1;
        }
        else {
            self.scale *= 0.9;
        }
        self.draw();
    }

    self.event_mouseover = function(e) {
        if(self.mfps > 0) {
            self.animate(self.mfps);
        }
    }

    self.event_mouseout = function(e) {
        if(self.mfps > 0) {
            self.animate(self.fps);
        }
    }


    self.queue_atoms = function() {
        for(var i = 0; i < self.atoms[self.index].length; i++) {
            var atom = self.atoms[self.index][i];
            var radius = elements[atom['symbol']]['radius'];
            var color = elements[atom['symbol']]['color'].slice(0);
            var position = atom['position'].slice(0);
            self.queue.push(["atom", position, radius, color, 0]);
        }

    }


    self.queue_transform = function() {
        for(var i = 0; i < self.queue.length; i++) {
            if(self.queue[i][0] == "atom") {
                self.queue[i][1][0] += self.translation[0];
                self.queue[i][1][1] += self.translation[1];
                self.queue[i][1][2] += self.translation[2];
                self.queue[i][1] = self.dot(self.rotation, self.queue[i][1]);
                self.queue[i][4] = self.queue[i][1][2];
            }
        }
    }


    self.queue_sort = function() {
        function sortby(a, b) {
            a = a[a.length - 1];
            b = b[b.length - 1];
            return a - b;
        }
        self.queue = self.queue.sort(sortby);
    }


    self.draw = function() {
        self.canvas.width = self.canvas.width;
        self.canvas.height = self.canvas.width;
        self.clear();
        self.queue = [];
        self.queue_atoms();
        self.queue_transform();
        self.queue_sort();
        self.draw_queue();
    }


    self.draw_queue = function() {
        for(var i = 0; i < self.queue.length; i++) {
            var q = self.queue[i];
            if(q[0] == "atom") {
                self.circle(Math.floor(q[1][0] * self.scale + self.canvas.width * 0.5), 
                            Math.floor(q[1][1] * self.scale + self.canvas.height * 0.5), 
                            Math.floor(q[2] * self.scale * 0.5 * self.radius), q[3]);
            }
        }
    }


    self.rgb = function(color) {
        var r = Math.floor(color[0] * 255.0);
        var g = Math.floor(color[1] * 255.0);
        var b = Math.floor(color[2] * 255.0);
        return "rgb(" + r + ", " + g + ", " + b + ")";
    }


    self.clear = function() {
        self.ctx.fillStyle = self.rgb(self.clearColor);
        self.ctx.fillRect(0, 0, self.canvas.width, self.canvas.height);
    }


    self.circle = function(x, y, radius, color) {
        self.ctx.lineWidth = 1.25;
        self.ctx.beginPath();
        color[0] = (color[0]) * 0.9;
        color[1] = (color[1]) * 0.9;
        color[2] = (color[2]) * 0.9;
        self.ctx.fillStyle = self.rgb(color);
        self.ctx.strokeStyle = "rgb(0, 0, 0)";
        self.ctx.arc(x, y, radius, 0, 2 * Math.PI, false);
        self.ctx.fill();
        self.ctx.stroke();
        self.ctx.beginPath();
        color[0] = (color[0] / 0.9 + 2.0) / 3.0;
        color[1] = (color[1] / 0.9 + 2.0) / 3.0;
        color[2] = (color[2] / 0.9 + 2.0) / 3.0;
        self.ctx.fillStyle = self.rgb(color);
        self.ctx.arc(x - radius / 4, y - radius / 4, radius / 4, 0, 2 * Math.PI, false);
        self.ctx.fill();
    }


    self.rotx = function(theta) {
        var ct = Math.cos(theta);
        var st = Math.sin(theta);
        var m = [[1, 0, 0], [0, ct, -st], [0, st, ct]];
        self.rotation = self.dot(m, self.rotation);
    }
    self.roty = function(theta) {
        var ct = Math.cos(theta);
        var st = Math.sin(theta);
        var m = [[ct, 0, st], [0, 1, 0], [-st, 0, ct]];
        self.rotation = self.dot(m, self.rotation);
    }
    self.rotz = function(theta) {
        var ct = Math.cos(theta);
        var st = Math.sin(theta);
        var m = [[ct, -st, 0], [st, ct, 0], [0, 0, 1]];
        self.rotation = self.dot(m, self.rotation);
    }


    self.dot = function(a, b) {
        if(b[0].length == 3) {
            return [[a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0], 
                     a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1], 
                     a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2]],
                    [a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0], 
                     a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1], 
                     a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2]],
                    [a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0], 
                     a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1], 
                     a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2]]];
        }
        else if(b[0].length == undefined) {
            return [a[0][0]*b[0]+a[0][1]*b[1]+a[0][2]*b[2], 
                    a[1][0]*b[0]+a[1][1]*b[1]+a[1][2]*b[2], 
                    a[2][0]*b[0]+a[2][1]*b[1]+a[2][2]*b[2]];
        }
    }

    self.setAtoms = function(a) {
        self.atoms = a.slice();
        atoms = self.atoms[0];
        var minx=1e200, maxx=-1e200, miny=1e200, maxy=-1e200, minz=1e200, maxz=-1e200;
        for(var i = 0; i < atoms.length; i++) {
            minx = Math.min(minx, atoms[i]['position'][0]);
            maxx = Math.max(maxx, atoms[i]['position'][0]);
            miny = Math.min(miny, atoms[i]['position'][1]);
            maxy = Math.max(maxy, atoms[i]['position'][1]);
            minz = Math.min(minz, atoms[i]['position'][2]);
            maxz = Math.max(maxz, atoms[i]['position'][2]);
        }
        self.translation[0] = -minx - ((maxx - minx) * 0.5);
        self.translation[1] = -miny - ((maxy - miny) * 0.5);
        self.translation[2] = -minz - ((maxz - minz) * 0.5);
        var xscale = self.canvas.width / ((maxx - minx) + 2*self.radius);
        var yscale = self.canvas.height / ((maxy - miny) + 2*self.radius);
        self.scale = Math.min(xscale, yscale);
        self.draw();
    }


    self.load = function(url, filetype, radius) {
        if(filetype == 'json') {
            self.loadJSON(url);
        }
        else if(filetype == 'xyz') {
            self.loadXYZ(url);
        }
        if(radius != undefined) {
            self.radius = radius;
        }
        self.animate(self.fps);
    }

    self.loadJSON = function(url) {
        jx.load(url, function(data) {
            self.setAtoms(JSON.parse(data));
        });
    }
    
    self.trim = function(str) {
        return str.replace(/^\s\s*/, '').replace(/\s\s*$/, '');
    }

    self.loadXYZ = function(url) {
        jx.load(url, function(data) {
            var lines = data.split('\n');
            var natoms = parseInt(lines[0]);
            var nframes = Math.floor(lines.length/(natoms+2));
            var trajectory = []
            for(var i = 0; i < nframes; i++) {
                var atoms = [];
                for(var j = 0; j < natoms; j++) {
                    var line = self.trim(lines[i*(natoms+2)+j+2]).split(/\s+/);
                    var atom = {};
                    atom.symbol = line[0];
                    atom.position = [parseFloat(line[1]), parseFloat(line[2]), parseFloat(line[3])];
                    atoms.push(atom);
                }
                trajectory.push(atoms);
            }
            self.setAtoms(trajectory);
        });
        
    }

    self.nextFrame = function(){
        self.index += 1;
        if(self.index >= self.atoms.length) {
            self.index = 0;
        }
        self.draw();
    }

    self.animate = function(fps) {
        clearInterval(self.animator);
        if (fps <= 0) {
            clearInterval(self.animator);
            return;
        }
        self.animator = setInterval(self.nextFrame, 1000/fps);
    }


    self.canvas = canvas;
    self.scale = 64.0;
    self.atoms = [];
    self.queue = [];
    self.rotation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    self.translation = [0, 0, 0];
    self.button1 = false;
    self.lastmousex = undefined;
    self.lastmousey = undefined;


    if(self.canvas == null) {
        return;
    }
    

    self.ctx = self.canvas.getContext("2d");
    window.addEventListener('mousemove', self.event_mousemove);
    window.addEventListener('mouseup', self.event_mouseup);
    self.canvas.onmouseover = self.event_mouseover;
    self.canvas.onmouseout = self.event_mouseout;
    self.canvas.onmousedown = self.event_mousedown;
    self.canvas.onmousewheel = self.event_mousewheel;


}

function _xyz_populate() {
    window._xyzs = [];
    var canvases = document.getElementsByClassName('xyz');
    for(var i = 0; i < canvases.length; i++) {
        var canvas = canvases[i];
        var newXYZ = new xyz(canvas);
        window._xyzs.push(newXYZ);
        newXYZ.mfps = canvas.getAttribute('mfps') ? parseInt(canvas.getAttribute('mfps')) : 0;
        newXYZ.fps = canvas.getAttribute('fps') ? parseInt(canvas.getAttribute('fps')) : 0;
        newXYZ.load(canvas.getAttribute('url'), 
                    canvas.getAttribute('filetype'), 
                    canvas.getAttribute('radius'));
    }
}

document.addEventListener("DOMContentLoaded", _xyz_populate);


var elements = {};
elements[  0] = elements[ 'Xx'] = {'symbol':  'Xx', 'name':       'unknown', 'mass':   1.00000000, 'radius':  1.0000, 'color': [1.000, 0.078, 0.576], 'number': 0};
elements[  1] = elements[  'H'] = {'symbol':   'H', 'name':      'hydrogen', 'mass':   1.00794000, 'radius':  0.3100, 'color': [1.000, 1.000, 1.000], 'number': 1};
elements[  2] = elements[ 'He'] = {'symbol':  'He', 'name':        'helium', 'mass':   4.00260200, 'radius':  0.2800, 'color': [0.851, 1.000, 1.000], 'number': 2};
elements[  3] = elements[ 'Li'] = {'symbol':  'Li', 'name':       'lithium', 'mass':   6.94100000, 'radius':  1.2800, 'color': [0.800, 0.502, 1.000], 'number': 3};
elements[  4] = elements[ 'Be'] = {'symbol':  'Be', 'name':     'beryllium', 'mass':   9.01218200, 'radius':  0.9600, 'color': [0.761, 1.000, 0.000], 'number': 4};
elements[  5] = elements[  'B'] = {'symbol':   'B', 'name':         'boron', 'mass':  10.81100000, 'radius':  0.8400, 'color': [1.000, 0.710, 0.710], 'number': 5};
elements[  6] = elements[  'C'] = {'symbol':   'C', 'name':        'carbon', 'mass':  12.01070000, 'radius':  0.7300, 'color': [0.565, 0.565, 0.565], 'number': 6};
elements[  7] = elements[  'N'] = {'symbol':   'N', 'name':      'nitrogen', 'mass':  14.00670000, 'radius':  0.7100, 'color': [0.188, 0.314, 0.973], 'number': 7};
elements[  8] = elements[  'O'] = {'symbol':   'O', 'name':        'oxygen', 'mass':  15.99940000, 'radius':  0.6600, 'color': [1.000, 0.051, 0.051], 'number': 8};
elements[  9] = elements[  'F'] = {'symbol':   'F', 'name':      'fluorine', 'mass':  18.99840320, 'radius':  0.5700, 'color': [0.565, 0.878, 0.314], 'number': 9};
elements[ 10] = elements[ 'Ne'] = {'symbol':  'Ne', 'name':          'neon', 'mass':  20.17970000, 'radius':  0.5800, 'color': [0.702, 0.890, 0.961], 'number': 10};
elements[ 11] = elements[ 'Na'] = {'symbol':  'Na', 'name':        'sodium', 'mass':  22.98976928, 'radius':  1.6600, 'color': [0.671, 0.361, 0.949], 'number': 11};
elements[ 12] = elements[ 'Mg'] = {'symbol':  'Mg', 'name':     'magnesium', 'mass':  24.30500000, 'radius':  1.4100, 'color': [0.541, 1.000, 0.000], 'number': 12};
elements[ 13] = elements[ 'Al'] = {'symbol':  'Al', 'name':      'aluminum', 'mass':  26.98153860, 'radius':  1.2100, 'color': [0.749, 0.651, 0.651], 'number': 13};
elements[ 14] = elements[ 'Si'] = {'symbol':  'Si', 'name':       'silicon', 'mass':  28.08550000, 'radius':  1.1100, 'color': [0.941, 0.784, 0.627], 'number': 14};
elements[ 15] = elements[  'P'] = {'symbol':   'P', 'name':    'phosphorus', 'mass':  30.97376200, 'radius':  1.0700, 'color': [1.000, 0.502, 0.000], 'number': 15};
elements[ 16] = elements[  'S'] = {'symbol':   'S', 'name':        'sulfur', 'mass':  32.06500000, 'radius':  1.0500, 'color': [1.000, 1.000, 0.188], 'number': 16};
elements[ 17] = elements[ 'Cl'] = {'symbol':  'Cl', 'name':      'chlorine', 'mass':  35.45300000, 'radius':  1.0200, 'color': [0.122, 0.941, 0.122], 'number': 17};
elements[ 18] = elements[ 'Ar'] = {'symbol':  'Ar', 'name':         'argon', 'mass':  39.94800000, 'radius':  1.0600, 'color': [0.502, 0.820, 0.890], 'number': 18};
elements[ 19] = elements[  'K'] = {'symbol':   'K', 'name':     'potassium', 'mass':  39.09830000, 'radius':  2.0300, 'color': [0.561, 0.251, 0.831], 'number': 19};
elements[ 20] = elements[ 'Ca'] = {'symbol':  'Ca', 'name':       'calcium', 'mass':  40.07800000, 'radius':  1.7600, 'color': [0.239, 1.000, 0.000], 'number': 20};
elements[ 21] = elements[ 'Sc'] = {'symbol':  'Sc', 'name':      'scandium', 'mass':  44.95591200, 'radius':  1.7000, 'color': [0.902, 0.902, 0.902], 'number': 21};
elements[ 22] = elements[ 'Ti'] = {'symbol':  'Ti', 'name':      'titanium', 'mass':  47.86700000, 'radius':  1.6000, 'color': [0.749, 0.761, 0.780], 'number': 22};
elements[ 23] = elements[  'V'] = {'symbol':   'V', 'name':      'vanadium', 'mass':  50.94150000, 'radius':  1.5300, 'color': [0.651, 0.651, 0.671], 'number': 23};
elements[ 24] = elements[ 'Cr'] = {'symbol':  'Cr', 'name':      'chromium', 'mass':  51.99610000, 'radius':  1.3900, 'color': [0.541, 0.600, 0.780], 'number': 24};
elements[ 25] = elements[ 'Mn'] = {'symbol':  'Mn', 'name':     'manganese', 'mass':  54.93804500, 'radius':  1.3900, 'color': [0.611, 0.478, 0.780], 'number': 25};
elements[ 26] = elements[ 'Fe'] = {'symbol':  'Fe', 'name':          'iron', 'mass':  55.84500000, 'radius':  1.3200, 'color': [0.878, 0.400, 0.200], 'number': 26};
elements[ 27] = elements[ 'Co'] = {'symbol':  'Co', 'name':        'cobalt', 'mass':  58.69340000, 'radius':  1.2600, 'color': [0.941, 0.565, 0.627], 'number': 27};
elements[ 28] = elements[ 'Ni'] = {'symbol':  'Ni', 'name':        'nickel', 'mass':  58.93319500, 'radius':  1.2400, 'color': [0.314, 0.816, 0.314], 'number': 28};
elements[ 29] = elements[ 'Cu'] = {'symbol':  'Cu', 'name':        'copper', 'mass':  63.54600000, 'radius':  1.3200, 'color': [0.784, 0.502, 0.200], 'number': 29};
elements[ 30] = elements[ 'Zn'] = {'symbol':  'Zn', 'name':          'zinc', 'mass':  65.38000000, 'radius':  1.2200, 'color': [0.490, 0.502, 0.690], 'number': 30};
elements[ 31] = elements[ 'Ga'] = {'symbol':  'Ga', 'name':       'gallium', 'mass':  69.72300000, 'radius':  1.2200, 'color': [0.761, 0.561, 0.561], 'number': 31};
elements[ 32] = elements[ 'Ge'] = {'symbol':  'Ge', 'name':     'germanium', 'mass':  72.64000000, 'radius':  1.2000, 'color': [0.400, 0.561, 0.561], 'number': 32};
elements[ 33] = elements[ 'As'] = {'symbol':  'As', 'name':       'arsenic', 'mass':  74.92160000, 'radius':  1.1900, 'color': [0.741, 0.502, 0.890], 'number': 33};
elements[ 34] = elements[ 'Se'] = {'symbol':  'Se', 'name':      'selenium', 'mass':  78.96000000, 'radius':  1.2000, 'color': [1.000, 0.631, 0.000], 'number': 34};
elements[ 35] = elements[ 'Br'] = {'symbol':  'Br', 'name':       'bromine', 'mass':  79.90400000, 'radius':  1.2000, 'color': [0.651, 0.161, 0.161], 'number': 35};
elements[ 36] = elements[ 'Kr'] = {'symbol':  'Kr', 'name':       'krypton', 'mass':  83.79800000, 'radius':  1.1600, 'color': [0.361, 0.722, 0.820], 'number': 36};
elements[ 37] = elements[ 'Rb'] = {'symbol':  'Rb', 'name':      'rubidium', 'mass':  85.46780000, 'radius':  2.2000, 'color': [0.439, 0.180, 0.690], 'number': 37};
elements[ 38] = elements[ 'Sr'] = {'symbol':  'Sr', 'name':     'strontium', 'mass':  87.62000000, 'radius':  1.9500, 'color': [0.000, 1.000, 0.000], 'number': 38};
elements[ 39] = elements[  'Y'] = {'symbol':   'Y', 'name':       'yttrium', 'mass':  88.90585000, 'radius':  1.9000, 'color': [0.580, 1.000, 1.000], 'number': 39};
elements[ 40] = elements[ 'Zr'] = {'symbol':  'Zr', 'name':     'zirconium', 'mass':  91.22400000, 'radius':  1.7500, 'color': [0.580, 0.878, 0.878], 'number': 40};
elements[ 41] = elements[ 'Nb'] = {'symbol':  'Nb', 'name':       'niobium', 'mass':  92.90638000, 'radius':  1.6400, 'color': [0.451, 0.761, 0.788], 'number': 41};
elements[ 42] = elements[ 'Mo'] = {'symbol':  'Mo', 'name':    'molybdenum', 'mass':  95.96000000, 'radius':  1.5400, 'color': [0.329, 0.710, 0.710], 'number': 42};
elements[ 43] = elements[ 'Tc'] = {'symbol':  'Tc', 'name':    'technetium', 'mass':  98.00000000, 'radius':  1.4700, 'color': [0.231, 0.620, 0.620], 'number': 43};
elements[ 44] = elements[ 'Ru'] = {'symbol':  'Ru', 'name':     'ruthenium', 'mass': 101.07000000, 'radius':  1.4600, 'color': [0.141, 0.561, 0.561], 'number': 44};
elements[ 45] = elements[ 'Rh'] = {'symbol':  'Rh', 'name':       'rhodium', 'mass': 102.90550000, 'radius':  1.4200, 'color': [0.039, 0.490, 0.549], 'number': 45};
elements[ 46] = elements[ 'Pd'] = {'symbol':  'Pd', 'name':     'palladium', 'mass': 106.42000000, 'radius':  1.3900, 'color': [0.000, 0.412, 0.522], 'number': 46};
elements[ 47] = elements[ 'Ag'] = {'symbol':  'Ag', 'name':        'silver', 'mass': 107.86820000, 'radius':  1.4500, 'color': [0.753, 0.753, 0.753], 'number': 47};
elements[ 48] = elements[ 'Cd'] = {'symbol':  'Cd', 'name':       'cadmium', 'mass': 112.41100000, 'radius':  1.4400, 'color': [1.000, 0.851, 0.561], 'number': 48};
elements[ 49] = elements[ 'In'] = {'symbol':  'In', 'name':        'indium', 'mass': 114.81800000, 'radius':  1.4200, 'color': [0.651, 0.459, 0.451], 'number': 49};
elements[ 50] = elements[ 'Sn'] = {'symbol':  'Sn', 'name':           'tin', 'mass': 118.71000000, 'radius':  1.3900, 'color': [0.400, 0.502, 0.502], 'number': 50};
elements[ 51] = elements[ 'Sb'] = {'symbol':  'Sb', 'name':      'antimony', 'mass': 121.76000000, 'radius':  1.3900, 'color': [0.620, 0.388, 0.710], 'number': 51};
elements[ 52] = elements[ 'Te'] = {'symbol':  'Te', 'name':     'tellurium', 'mass': 127.60000000, 'radius':  1.3800, 'color': [0.831, 0.478, 0.000], 'number': 52};
elements[ 53] = elements[  'I'] = {'symbol':   'I', 'name':        'iodine', 'mass': 126.90470000, 'radius':  1.3900, 'color': [0.580, 0.000, 0.580], 'number': 53};
elements[ 54] = elements[ 'Xe'] = {'symbol':  'Xe', 'name':         'xenon', 'mass': 131.29300000, 'radius':  1.4000, 'color': [0.259, 0.620, 0.690], 'number': 54};
elements[ 55] = elements[ 'Cs'] = {'symbol':  'Cs', 'name':        'cesium', 'mass': 132.90545190, 'radius':  2.4400, 'color': [0.341, 0.090, 0.561], 'number': 55};
elements[ 56] = elements[ 'Ba'] = {'symbol':  'Ba', 'name':        'barium', 'mass': 137.32700000, 'radius':  2.1500, 'color': [0.000, 0.788, 0.000], 'number': 56};
elements[ 57] = elements[ 'La'] = {'symbol':  'La', 'name':     'lanthanum', 'mass': 138.90547000, 'radius':  2.0700, 'color': [0.439, 0.831, 1.000], 'number': 57};
elements[ 58] = elements[ 'Ce'] = {'symbol':  'Ce', 'name':        'cerium', 'mass': 140.11600000, 'radius':  2.0400, 'color': [1.000, 1.000, 0.780], 'number': 58};
elements[ 59] = elements[ 'Pr'] = {'symbol':  'Pr', 'name':  'praseodymium', 'mass': 140.90765000, 'radius':  2.0300, 'color': [0.851, 1.000, 0.780], 'number': 59};
elements[ 60] = elements[ 'Nd'] = {'symbol':  'Nd', 'name':     'neodymium', 'mass': 144.24200000, 'radius':  2.0100, 'color': [0.780, 1.000, 0.780], 'number': 60};
elements[ 61] = elements[ 'Pm'] = {'symbol':  'Pm', 'name':    'promethium', 'mass': 145.00000000, 'radius':  1.9900, 'color': [0.639, 1.000, 0.780], 'number': 61};
elements[ 62] = elements[ 'Sm'] = {'symbol':  'Sm', 'name':      'samarium', 'mass': 150.36000000, 'radius':  1.9800, 'color': [0.561, 1.000, 0.780], 'number': 62};
elements[ 63] = elements[ 'Eu'] = {'symbol':  'Eu', 'name':      'europium', 'mass': 151.96400000, 'radius':  1.9800, 'color': [0.380, 1.000, 0.780], 'number': 63};
elements[ 64] = elements[ 'Gd'] = {'symbol':  'Gd', 'name':    'gadolinium', 'mass': 157.25000000, 'radius':  1.9600, 'color': [0.271, 1.000, 0.780], 'number': 64};
elements[ 65] = elements[ 'Tb'] = {'symbol':  'Tb', 'name':       'terbium', 'mass': 158.92535000, 'radius':  1.9400, 'color': [0.189, 1.000, 0.780], 'number': 65};
elements[ 66] = elements[ 'Dy'] = {'symbol':  'Dy', 'name':    'dysprosium', 'mass': 162.50000000, 'radius':  1.9200, 'color': [0.122, 1.000, 0.780], 'number': 66};
elements[ 67] = elements[ 'Ho'] = {'symbol':  'Ho', 'name':       'holmium', 'mass': 164.93032000, 'radius':  1.9200, 'color': [0.000, 1.000, 0.612], 'number': 67};
elements[ 68] = elements[ 'Er'] = {'symbol':  'Er', 'name':        'erbium', 'mass': 167.25900000, 'radius':  1.8900, 'color': [0.000, 0.902, 0.459], 'number': 68};
elements[ 69] = elements[ 'Tm'] = {'symbol':  'Tm', 'name':       'thulium', 'mass': 168.93421000, 'radius':  1.9000, 'color': [0.000, 0.831, 0.322], 'number': 69};
elements[ 70] = elements[ 'Yb'] = {'symbol':  'Yb', 'name':     'ytterbium', 'mass': 173.05400000, 'radius':  1.8700, 'color': [0.000, 0.749, 0.220], 'number': 70};
elements[ 71] = elements[ 'Lu'] = {'symbol':  'Lu', 'name':      'lutetium', 'mass': 174.96680000, 'radius':  1.8700, 'color': [0.000, 0.671, 0.141], 'number': 71};
elements[ 72] = elements[ 'Hf'] = {'symbol':  'Hf', 'name':       'hafnium', 'mass': 178.49000000, 'radius':  1.7500, 'color': [0.302, 0.761, 1.000], 'number': 72};
elements[ 73] = elements[ 'Ta'] = {'symbol':  'Ta', 'name':      'tantalum', 'mass': 180.94788000, 'radius':  1.7000, 'color': [0.302, 0.651, 1.000], 'number': 73};
elements[ 74] = elements[  'W'] = {'symbol':   'W', 'name':      'tungsten', 'mass': 183.84000000, 'radius':  1.6200, 'color': [0.129, 0.580, 0.839], 'number': 74};
elements[ 75] = elements[ 'Re'] = {'symbol':  'Re', 'name':       'rhenium', 'mass': 186.20700000, 'radius':  1.5100, 'color': [0.149, 0.490, 0.671], 'number': 75};
elements[ 76] = elements[ 'Os'] = {'symbol':  'Os', 'name':        'osmium', 'mass': 190.23000000, 'radius':  1.4400, 'color': [0.149, 0.400, 0.588], 'number': 76};
elements[ 77] = elements[ 'Ir'] = {'symbol':  'Ir', 'name':       'iridium', 'mass': 192.21700000, 'radius':  1.4100, 'color': [0.090, 0.329, 0.529], 'number': 77};
elements[ 78] = elements[ 'Pt'] = {'symbol':  'Pt', 'name':      'platinum', 'mass': 195.08400000, 'radius':  1.3600, 'color': [0.816, 0.816, 0.878], 'number': 78};
elements[ 79] = elements[ 'Au'] = {'symbol':  'Au', 'name':          'gold', 'mass': 196.96656900, 'radius':  1.3600, 'color': [1.000, 0.820, 0.137], 'number': 79};
elements[ 80] = elements[ 'Hg'] = {'symbol':  'Hg', 'name':       'mercury', 'mass': 200.59000000, 'radius':  1.3200, 'color': [0.722, 0.722, 0.816], 'number': 80};
elements[ 81] = elements[ 'Tl'] = {'symbol':  'Tl', 'name':      'thallium', 'mass': 204.38330000, 'radius':  1.4500, 'color': [0.651, 0.329, 0.302], 'number': 81};
elements[ 82] = elements[ 'Pb'] = {'symbol':  'Pb', 'name':          'lead', 'mass': 207.20000000, 'radius':  1.4600, 'color': [0.341, 0.349, 0.380], 'number': 82};
elements[ 83] = elements[ 'Bi'] = {'symbol':  'Bi', 'name':       'bismuth', 'mass': 208.98040000, 'radius':  1.4800, 'color': [0.620, 0.310, 0.710], 'number': 83};
elements[ 84] = elements[ 'Po'] = {'symbol':  'Po', 'name':      'polonium', 'mass': 210.00000000, 'radius':  1.4000, 'color': [0.671, 0.361, 0.000], 'number': 84};
elements[ 85] = elements[ 'At'] = {'symbol':  'At', 'name':      'astatine', 'mass': 210.00000000, 'radius':  1.5000, 'color': [0.459, 0.310, 0.271], 'number': 85};
elements[ 86] = elements[ 'Rn'] = {'symbol':  'Rn', 'name':         'radon', 'mass': 220.00000000, 'radius':  1.5000, 'color': [0.259, 0.510, 0.588], 'number': 86};
elements[ 87] = elements[ 'Fr'] = {'symbol':  'Fr', 'name':      'francium', 'mass': 223.00000000, 'radius':  2.6000, 'color': [0.259, 0.000, 0.400], 'number': 87};
elements[ 88] = elements[ 'Ra'] = {'symbol':  'Ra', 'name':        'radium', 'mass': 226.00000000, 'radius':  2.2100, 'color': [0.000, 0.490, 0.000], 'number': 88};
elements[ 89] = elements[ 'Ac'] = {'symbol':  'Ac', 'name':      'actinium', 'mass': 227.00000000, 'radius':  2.1500, 'color': [0.439, 0.671, 0.980], 'number': 89};
elements[ 90] = elements[ 'Th'] = {'symbol':  'Th', 'name':       'thorium', 'mass': 231.03588000, 'radius':  2.0600, 'color': [0.000, 0.729, 1.000], 'number': 90};
elements[ 91] = elements[ 'Pa'] = {'symbol':  'Pa', 'name':  'protactinium', 'mass': 232.03806000, 'radius':  2.0000, 'color': [0.000, 0.631, 1.000], 'number': 91};
elements[ 92] = elements[  'U'] = {'symbol':   'U', 'name':       'uranium', 'mass': 237.00000000, 'radius':  1.9600, 'color': [0.000, 0.561, 1.000], 'number': 92};
elements[ 93] = elements[ 'Np'] = {'symbol':  'Np', 'name':     'neptunium', 'mass': 238.02891000, 'radius':  1.9000, 'color': [0.000, 0.502, 1.000], 'number': 93};
elements[ 94] = elements[ 'Pu'] = {'symbol':  'Pu', 'name':     'plutonium', 'mass': 243.00000000, 'radius':  1.8700, 'color': [0.000, 0.420, 1.000], 'number': 94};
elements[ 95] = elements[ 'Am'] = {'symbol':  'Am', 'name':     'americium', 'mass': 244.00000000, 'radius':  1.8000, 'color': [0.329, 0.361, 0.949], 'number': 95};
elements[ 96] = elements[ 'Cm'] = {'symbol':  'Cm', 'name':        'curium', 'mass': 247.00000000, 'radius':  1.6900, 'color': [0.471, 0.361, 0.890], 'number': 96};
elements[ 97] = elements[ 'Bk'] = {'symbol':  'Bk', 'name':     'berkelium', 'mass': 247.00000000, 'radius':  1.6600, 'color': [0.541, 0.310, 0.890], 'number': 97};
elements[ 98] = elements[ 'Cf'] = {'symbol':  'Cf', 'name':   'californium', 'mass': 251.00000000, 'radius':  1.6800, 'color': [0.631, 0.212, 0.831], 'number': 98};
elements[ 99] = elements[ 'Es'] = {'symbol':  'Es', 'name':   'einsteinium', 'mass': 252.00000000, 'radius':  1.6500, 'color': [0.702, 0.122, 0.831], 'number': 99};
elements[100] = elements[ 'Fm'] = {'symbol':  'Fm', 'name':       'fermium', 'mass': 257.00000000, 'radius':  1.6700, 'color': [0.702, 0.122, 0.729], 'number': 100};
elements[101] = elements[ 'Md'] = {'symbol':  'Md', 'name':   'mendelevium', 'mass': 258.00000000, 'radius':  1.7300, 'color': [0.702, 0.051, 0.651], 'number': 101};
elements[102] = elements[ 'No'] = {'symbol':  'No', 'name':      'nobelium', 'mass': 259.00000000, 'radius':  1.7600, 'color': [0.741, 0.051, 0.529], 'number': 102};
elements[103] = elements[ 'Lr'] = {'symbol':  'Lr', 'name':    'lawrencium', 'mass': 262.00000000, 'radius':  1.6100, 'color': [0.780, 0.000, 0.400], 'number': 103};
elements[104] = elements[ 'Rf'] = {'symbol':  'Rf', 'name': 'rutherfordium', 'mass': 261.00000000, 'radius':  1.5700, 'color': [0.800, 0.000, 0.349], 'number': 104};
elements[105] = elements[ 'Db'] = {'symbol':  'Db', 'name':       'dubnium', 'mass': 262.00000000, 'radius':  1.4900, 'color': [0.820, 0.000, 0.310], 'number': 105};
elements[106] = elements[ 'Sg'] = {'symbol':  'Sg', 'name':    'seaborgium', 'mass': 266.00000000, 'radius':  1.4300, 'color': [0.851, 0.000, 0.271], 'number': 106};
elements[107] = elements[ 'Bh'] = {'symbol':  'Bh', 'name':       'bohrium', 'mass': 264.00000000, 'radius':  1.4100, 'color': [0.878, 0.000, 0.220], 'number': 107};
elements[108] = elements[ 'Hs'] = {'symbol':  'Hs', 'name':       'hassium', 'mass': 277.00000000, 'radius':  1.3400, 'color': [0.902, 0.000, 0.180], 'number': 108};
elements[109] = elements[ 'Mt'] = {'symbol':  'Mt', 'name':    'meitnerium', 'mass': 268.00000000, 'radius':  1.2900, 'color': [0.922, 0.000, 0.149], 'number': 109};
elements[110] = elements[ 'Ds'] = {'symbol':  'Ds', 'name':            'Ds', 'mass': 271.00000000, 'radius':  1.2800, 'color': [0.922, 0.000, 0.149], 'number': 110};
elements[111] = elements['Uuu'] = {'symbol': 'Uuu', 'name':           'Uuu', 'mass': 272.00000000, 'radius':  1.2100, 'color': [0.922, 0.000, 0.149], 'number': 111};
elements[112] = elements['Uub'] = {'symbol': 'Uub', 'name':           'Uub', 'mass': 285.00000000, 'radius':  1.2200, 'color': [0.922, 0.000, 0.149], 'number': 112};
elements[113] = elements['Uut'] = {'symbol': 'Uut', 'name':           'Uut', 'mass': 284.00000000, 'radius':  1.3600, 'color': [0.922, 0.000, 0.149], 'number': 113};
elements[114] = elements['Uuq'] = {'symbol': 'Uuq', 'name':           'Uuq', 'mass': 289.00000000, 'radius':  1.4300, 'color': [0.922, 0.000, 0.149], 'number': 114};
elements[115] = elements['Uup'] = {'symbol': 'Uup', 'name':           'Uup', 'mass': 288.00000000, 'radius':  1.6200, 'color': [0.922, 0.000, 0.149], 'number': 115};
elements[116] = elements['Uuh'] = {'symbol': 'Uuh', 'name':           'Uuh', 'mass': 292.00000000, 'radius':  1.7500, 'color': [0.922, 0.000, 0.149], 'number': 116};
elements[117] = elements['Uus'] = {'symbol': 'Uus', 'name':           'Uus', 'mass': 294.00000000, 'radius':  1.6500, 'color': [0.922, 0.000, 0.149], 'number': 117};
elements[118] = elements['Uuo'] = {'symbol': 'Uuo', 'name':           'Uuo', 'mass': 296.00000000, 'radius':  1.5700, 'color': [0.922, 0.000, 0.149], 'number': 118};




















