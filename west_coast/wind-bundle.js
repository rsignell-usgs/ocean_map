/**
 * Wind map code (c) 2012
 * Fernanda Viegas & Martin Wattenberg
 */

 /**
 * Simple representation of 2D vector.
 */

var Vector = function(x, y) {
	this.x = x;
	this.y = y;
}


Vector.polar = function(r, theta) {
	return new Vector(r * Math.cos(theta), r * Math.sin(theta));
};


Vector.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y);
};


Vector.prototype.copy = function(){
  return new Vector(this.x, this.y);
};


Vector.prototype.setLength = function(length) {
	var current = this.length();
	if (current) {
		var scale = length / current;
		this.x *= scale;
		this.y *= scale;
	}
	return this;
};


Vector.prototype.setAngle = function(theta) {
  var r = length();
  this.x = r * Math.cos(theta);
  this.y = r * Math.sin(theta);
  return this;
};


Vector.prototype.getAngle = function() {
  return Math.atan2(this.y, this.x);
};


Vector.prototype.d = function(v) {
		var dx = v.x - this.x;
		var dy = v.y - this.y;
		return Math.sqrt(dx * dx + dy * dy);
};/**
 * Identity projection.
 */
var IDProjection = {
	project: function(x, y, opt_v) {
		var v = opt_v || new Vector();
		v.x = x;
		v.y = y;
	  return v;
  },
	invert: function(x, y, opt_v) {
		var v = opt_v || new Vector();
		v.x = x;
		v.y = y;
	  return v;
  }
};

/**
 * Albers equal-area projection.
 * Constant param values after d3 (Bostock, Carden).
 */
var Albers = function() {
  function radians(degrees) {
		return Math.PI * degrees / 180;
  }

  var phi1 = radians(25.0);
  var phi2 = radians(52.0);
  var n = .5 * (phi1 + phi2);
	var C = Math.cos(phi1) * Math.cos(phi1) + 2 * n * Math.sin(phi1);
	var phi0 = radians(37.0);
	var lambda0 = radians(-125.0);
	var rho0 = Math.sqrt(C - 2 * n * Math.sin(phi0)) / n;

  return {
		project: function(lon, lat, opt_result) {
			lon = radians(lon);
		  lat = radians(lat);
		  var theta = n * (lon - lambda0);
		  var rho = Math.sqrt(C - 2 * n * Math.sin(lat)) / n;
		  var x = rho * Math.sin(theta);
		  var y = rho0 - rho * Math.cos(theta);
			if (opt_result) {
		    opt_result.x = x;
		    opt_result.y = y;
		    return opt_result;
	    }
		  return new Vector(x, y);
		},
		invert: function(x, y) {
			var rho2 = x * x + (rho0 - y) * (rho0 - y);
			var theta = Math.atan(x / (rho0 - y));
			var lon = lambda0 + theta / n;
			var lat = Math.asin((C / n - rho2 * n) / 2);
			return new Vector(lon * 180 / Math.PI, lat * 180 / Math.PI);
		}
	};
}();


var ScaledAlbers = function(scale, offsetX, offsetY, longMin, latMin) {
	this.scale = scale;
	this.offsetX = offsetX;
	this.offsetY = offsetY;
	this.longMin = longMin;
	this.latMin = latMin;
  this.swCorner = Albers.project(longMin, latMin);
};

ScaledAlbers.temp = new Vector(0, 0);

ScaledAlbers.prototype.project = function(lon, lat, opt_result) {
  var proj = Albers.project(lon, lat, ScaledAlbers.temp);
  var a = proj.x;
	var b = proj.y;
	var x = this.scale * (a - this.swCorner.x) + this.offsetX;
	var y = -this.scale * (b - this.swCorner.y) + this.offsetY;
	if (opt_result) {
		opt_result.x = x;
		opt_result.y = y;
		return opt_result;
	}
	return new Vector(x, y);
};

ScaledAlbers.prototype.invert = function(x, y) {
	var a = (x - this.offsetX) / this.scale + this.swCorner.x;
	var b = (y - this.offsetY) / -this.scale + this.swCorner.y;
	return Albers.invert(a, b);
};

/**
 * Represents a vector field based on an array of data,
 * with specified grid coordinates, using bilinear interpolation
 * for values that don't lie on grid points.
 */

/**
 * 
 * @param field 2D array of Vectors
 * 
 * next params are corners of region.
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 */
var VectorField = function(field, x0, y0, x1, y1) {
	this.x0 = x0;
	this.x1 = x1;
	this.y0 = y0;
	this.y1 = y1;
	this.field = field;
	this.w = field.length;
	this.h = field[0].length;
	this.maxLength = 0;
	var mx = 0;
	var my = 0;
	for (var i = 0; i < this.w; i++) {
	  for (var j = 0; j < this.h; j++) {
			if (field[i][j].length() > this.maxLength) {
				mx = i;
				my = j;
			}
			this.maxLength = Math.max(this.maxLength, field[i][j].length());
		}
	}
	mx = (mx / this.w) * (x1 - x0) + x0;
	my = (my / this.h) * (y1 - y0) + y0;
};

/**
 * Reads data from raw object in form:
 * {
 *   x0: -126.292942,
 *   y0: 23.525552,
 *   x1: -66.922962,
 *   y1: 49.397231,
 *   gridWidth: 501.0,
 *   gridHeight: 219.0,
 *   field: [
 *     0,0,
 *     0,0,
 *     ... (list of vectors)
 *   ]
 * }
 *
 * If the correctForSphere flag is set, we correct for the
 * distortions introduced by an equirectangular projection.
 */
VectorField.read = function(data, correctForSphere) {
	var field = [];
	var w = data.gridWidth;
	var h = data.gridHeight;
	var n = 2 * w * h;
	var i = 0;
	// OK, "total" and "weight"
	// are kludges that you should totally ignore,
	// unless you are interested in the average
	// vector length on vector field over lat/lon domain.
	var total = 0;
	var weight = 0;
	for (var x = 0; x < w; x++) {
		field[x] = [];
		for (var y = 0; y < h; y++) {
			var vx = data.field[i++] * 10;
			var vy = data.field[i++] * 10;
			var v = new Vector(vx, vy);
			// Uncomment to test a constant field:
			// v = new Vector(10, 0);
			if (correctForSphere) {
				var ux = x / (w - 1);
				var uy = y / (h - 1);
				var lon = data.x0 * (1 - ux) + data.x1 * ux;
				var lat = data.y0 * (1 - uy) + data.y1 * uy;
				var m = Math.PI * lat / 180;
				var length = v.length();
				if (length) {
			    total += length * m;
			    weight += m;
		    }
				v.x /= Math.cos(m);
				v.setLength(length);
			}
			field[x][y] = v;
		}
	}
	var result = new VectorField(field, data.x0, data.y0, data.x1, data.y1);
  //window.console.log('total = ' + total);
	//window.console.log('weight = ' + weight);
  if (total && weight) {

	  result.averageLength = total / weight;
	}
	return result;
};
  
VectorField.prototype.inBounds = function(x, y) {
  return x >= this.x0 && x < this.x1 && y >= this.y0 && y < this.y1;
};


VectorField.prototype.bilinear = function(coord, a, b) {
  var na = Math.floor(a);
  var nb = Math.floor(b);
  var ma = Math.ceil(a);
  var mb = Math.ceil(b);
  var fa = a - na;
  var fb = b - nb;

  return this.field[na][nb][coord] * (1 - fa) * (1 - fb) +
  	     this.field[ma][nb][coord] * fa * (1 - fb) +
  	     this.field[na][mb][coord] * (1 - fa) * fb +
  	     this.field[ma][mb][coord] * fa * fb;
};


VectorField.prototype.getValue = function(x, y, opt_result) {
	var a = (this.w - 1 - 1e-6) * (x - this.x0) / (this.x1 - this.x0);
	var b = (this.h - 1 - 1e-6) * (y - this.y0) / (this.y1 - this.y0);
	var vx = this.bilinear('x', a, b);
	var vy = this.bilinear('y', a, b);
	if (opt_result) {
		opt_result.x = vx;
		opt_result.y = vy;
		return opt_result;
	}
	return new Vector(vx, vy);
};


VectorField.prototype.vectValue = function(vector) {
	return this.getValue(vector.x, vector.y);
};


VectorField.constant = function(dx, dy, x0, y0, x1, y1) {
	var field = new VectorField([[]], x0, y0, x1, y1);
	field.maxLength = Math.sqrt(dx * dx + dy * dy);
	field.getValue = function() {
		return new Vector(dx, dy);
	}
	return field;
}
/**
 * Listens to mouse events on an element, tracks zooming and panning,
 * informs other components of what's going on.
 */
var Animator = function(element, opt_animFunc, opt_unzoomButton) {
 	this.element = element;
	this.mouseIsDown = false;
	this.mouseX = -1;
	this.mouseY = -1;
	this.animating = true;
	this.state = 'animate';
	this.listeners = [];
	this.dx = 0;
	this.dy = 0;
	this.scale = 1;
	this.zoomProgress = 0;
	this.scaleTarget = 1;
	this.scaleStart = 1;
	this.animFunc = opt_animFunc;
	this.unzoomButton = opt_unzoomButton;
	
	if (element) {
		var self = this;
  	$(element).mousedown(function(e){
			self.mouseX = e.pageX - this.offsetLeft;
	    self.mouseY = e.pageY - this.offsetTop;
  		self.mousedown();
  	});
  	$(element).mouseup(function(e){
			self.mouseX = e.pageX - this.offsetLeft;
	    self.mouseY = e.pageY - this.offsetTop;
  		self.mouseup();
  	});
  	$(element).mousemove(function(e){
			self.mouseX = e.pageX - this.offsetLeft;
	    self.mouseY = e.pageY - this.offsetTop;
  		self.mousemove();
  	});
  }
};
 

Animator.prototype.mousedown = function() {
	this.state = 'mouse-down';
	this.notify('startMove');
	this.landingX = this.mouseX;
	this.landingY = this.mouseY;
	this.dxStart = this.dx;
	this.dyStart = this.dy;
	this.scaleStart = this.scale;
	this.mouseIsDown = true;
};


Animator.prototype.mousemove = function() {
	if (!this.mouseIsDown) {
		this.notify('hover');
		return;
	}
	var ddx = this.mouseX - this.landingX;
	var ddy = this.mouseY - this.landingY;
	var slip = Math.abs(ddx) + Math.abs(ddy);
	if (slip > 2 || this.state == 'pan') {
		this.state = 'pan';
		this.dx += ddx;
		this.dy += ddy;
		this.landingX = this.mouseX;
		this.landingY = this.mouseY;
		this.notify('move');
	}
}

Animator.prototype.mouseup = function() {
	this.mouseIsDown = false;
	if (this.state == 'pan') {
		this.state = 'animate';
		this.notify('endMove');
		return;
	}
	this.zoomClick(this.mouseX, this.mouseY);
};

 
Animator.prototype.add = function(listener) {
 	this.listeners.push(listener);
};


Animator.prototype.notify = function(message) {
	if (this.unzoomButton) {
		var diff = Math.abs(this.scale - 1) > .001 ||
		           Math.abs(this.dx) > .001 || Math.abs(this.dy > .001);
		this.unzoomButton.style.visibility = diff ? 'visible' : 'hidden';
	}
	if (this.animFunc && !this.animFunc()) {
		return;
	}
	for (var i = 0; i < this.listeners.length; i++) {
		var listener = this.listeners[i];
		if (listener[message]) {
			listener[message].call(listener, this);
		}
	}
};


Animator.prototype.unzoom = function() {
	this.zoom(0, 0, 1);
};


Animator.prototype.zoomClick = function(x, y) {
	var z = 1.7;
	var scale = 1.7 * this.scale;
	var dx = x - z * (x - this.dx);
	var dy = y - z * (y - this.dy);
	this.zoom(dx, dy, scale);
};

Animator.prototype.zoom = function(dx, dy, scale) {
	this.state = 'zoom';
  this.zoomProgress = 0;
  this.scaleStart = this.scale;
	this.scaleTarget = scale;
	this.dxTarget = dx;
	this.dyTarget = dy;
	this.dxStart = this.dx;
	this.dyStart = this.dy;
	this.notify('startMove');
};

Animator.prototype.relativeZoom = function() {
	return this.scale / this.scaleStart;
};


Animator.prototype.relativeDx = function() {
	return this.dx - this.dxStart;
}

Animator.prototype.relativeDy = function() {
	return this.dy - this.dyStart;
}

Animator.prototype.start = function(opt_millis) {
	var millis = opt_millis || 20;
	var self = this;
	function go() {
		var start = new Date();
		self.loop();
		var time = new Date() - start;
		setTimeout(go, Math.max(10, millis - time));
	}
	go();
};


Animator.prototype.loop = function() {
	if (this.state == 'mouse-down' || this.state == 'pan') {
		return;
	}
	if (this.state == 'animate') {
  	this.notify('animate');
		return;
  }
	if (this.state == 'zoom') {
  	this.zoomProgress = Math.min(1, this.zoomProgress + .07);
	  var u = (1 + Math.cos(Math.PI * this.zoomProgress)) / 2;
		function lerp(a, b) {
			return u * a + (1 - u) * b;
		}
	  this.scale = lerp(this.scaleStart, this.scaleTarget);
		this.dx = lerp(this.dxStart, this.dxTarget);
		this.dy = lerp(this.dyStart, this.dyTarget);
  	if (this.zoomProgress < 1) {
  		this.notify('move');
  	} else {
  		this.state = 'animate';
  		this.zoomCurrent = this.zoomTarget;
   		this.notify('endMove');
  	}
  }
};
 
/**
 * Displays a geographic vector field using moving particles.
 * Positions in the field are drawn onscreen using the Alber
 * "Projection" file.
 */

var Particle = function(x, y, age) {
	this.x = x;
	this.y = y;
	this.oldX = -1;
	this.oldY = -1;
	this.age = age;
	this.rnd = Math.random();
}


/**
 * @param {HTMLCanvasElement} canvas
 * @param {number} scale The scale factor for the projection.
 * @param {number} offsetX
 * @param {number} offsetY
 * @param {number} longMin
 * @param {number} latMin
 * @param {VectorField} field
 * @param {number} numParticles
 */
var MotionDisplay = function(canvas, imageCanvas, field, numParticles, opt_projection) {
	this.canvas = canvas;
  this.projection = opt_projection || IDProjection;
  this.field = field;
	this.numParticles = numParticles;
	this.first = true;
	this.maxLength = field.maxLength;
	this.speedScale = 2.5;   // controls the speed of the streaks
	this.renderState = 'normal';
	this.imageCanvas = imageCanvas;
	this.x0 = this.field.x0;
	this.x1 = this.field.x1;
	this.y0 = this.field.y0;
	this.y1 = this.field.y1;
	this.makeNewParticles(null, true);
	this.colors = [];
	this.rgb = '40, 40, 40';
	this.background = 'rgb(' + this.rgb + ')';
	this.backgroundAlpha = 'rgba(' + this.rgb + ', .02)';
	this.outsideColor = '#fff';
	for (var i = 0; i < 256; i++) {
		this.colors[i] = 'rgb(' + i + ',' + i + ',' + i + ')';
	}
	if (this.projection) {
  	this.startOffsetX = this.projection.offsetX;
  	this.startOffsetY = this.projection.offsetY;
  	this.startScale = this.projection.scale;
  }
};


MotionDisplay.prototype.setAlpha = function(alpha) {
	this.backgroundAlpha = 'rgba(' + this.rgb + ', ' + alpha + ')';
};

MotionDisplay.prototype.makeNewParticles = function(animator) {
	this.particles = [];
	for (var i = 0; i < this.numParticles; i++) {
		this.particles.push(this.makeParticle(animator));
	}
};


MotionDisplay.prototype.makeParticle = function(animator) {
	var dx = animator ? animator.dx : 0;
	var dy = animator ? animator.dy : 0;
	var scale = animator ? animator.scale : 1;
	var safecount = 0;
	for (;;) {
		var a = Math.random();
		var b = Math.random();
		var x = a * this.x0 + (1 - a) * this.x1;
		var y = b * this.y0 + (1 - b) * this.y1;
		var v = this.field.getValue(x, y);
		if (this.field.maxLength == 0) {
			return new Particle(x, y, 1 + 40 * Math.random());
		}
		var m = v.length() / this.field.maxLength;
		// The random factor here is designed to ensure that
		// more particles are placed in slower areas; this makes the
		// overall distribution appear more even.
		if ((v.x || v.y) && (++safecount > 10 || Math.random() > m * .9)) {
			var proj = this.projection.project(x, y);
			var sx = proj.x * scale + dx;
			var sy = proj.y * scale + dy;
			if (++safecount > 10 || !(sx < 0 || sy < 0 || sx > this.canvas.width || sy > this.canvas.height)) {
	      return new Particle(x, y, 1 + 40 * Math.random());
      }	
		}
	}
};


MotionDisplay.prototype.startMove = function(animator) {
	// Save screen.
	this.imageCanvas.getContext('2d').drawImage(this.canvas, 0, 0);
};


MotionDisplay.prototype.endMove  = function(animator) {
	if (animator.scale < 1.1) {
		this.x0 = this.field.x0;
		this.x1 = this.field.x1;
		this.y0 = this.field.y0;
		this.y1 = this.field.y1;
	} else {
		// get new bounds for making new particles.
		var p = this.projection;
		var self = this;
		function invert(x, y) {
			x = (x - animator.dx) / animator.scale;
			y = (y - animator.dy) / animator.scale;
			return self.projection.invert(x, y);
		}
		var loc = invert(0, 0);
		var x0 = loc.x;
		var x1 = loc.x;
		var y0 = loc.y;
		var y1 = loc.y;
		function expand(x, y) {
			var v = invert(x, y);
			x0 = Math.min(v.x, x0);
			x1 = Math.max(v.x, x1);
			y0 = Math.min(v.y, y0);
			y1 = Math.max(v.y, y1);
		}
		// This calculation with "top" is designed to fix a bug
		// where we were missing particles at the top of the
		// screen with north winds. This is a short-term fix,
		// it's dependent on the particular projection and
		// region, and we should figure out a more general
		// solution soon.
		var top = -.2 * this.canvas.height;
		expand(top, this.canvas.height);
		expand(this.canvas.width, top);
		expand(this.canvas.width, this.canvas.height);
		this.x0 = Math.max(this.field.x0, x0);
		this.x1 = Math.min(this.field.x1, x1);
		this.y0 = Math.max(this.field.y0, y0);
		this.y1 = Math.min(this.field.y1, y1);
	}
	tick = 0;
	this.makeNewParticles(animator);
};


MotionDisplay.prototype.animate = function(animator) {
	this.moveThings(animator);
  this.draw(animator);
}


MotionDisplay.prototype.move = function(animator) {
	var w = this.canvas.width;
	var h = this.canvas.height;
	var g = this.canvas.getContext('2d');
	
	g.fillStyle = this.outsideColor;
	var dx = animator.dx;
	var dy = animator.dy;
	var scale = animator.scale;

	g.fillRect(0, 0, w, h);
	g.fillStyle = this.background;
  g.fillRect(dx, dy, w * scale, h * scale);
	var z = animator.relativeZoom();
	var dx = animator.dx - z * animator.dxStart;
	var dy = animator.dy - z * animator.dyStart;
	g.drawImage(this.imageCanvas, dx, dy, z * w, z * h);
};


MotionDisplay.prototype.moveThings = function(animator) {
	var speed = .01 * this.speedScale / animator.scale;
	for (var i = 0; i < this.particles.length; i++) {
		var p = this.particles[i];
		if (p.age > 0 && this.field.inBounds(p.x, p.y)) {
		  var a = this.field.getValue(p.x, p.y);
			p.x += speed * a.x;
			p.y += speed * a.y;
			p.age--;
		} else {
			this.particles[i] = this.makeParticle(animator);
		}
	}
};


MotionDisplay.prototype.draw = function(animator) {
	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	if (this.first) {
		g.fillStyle =  this.background;
		this.first = false;
	} else {
		g.fillStyle = this.backgroundAlpha;
	}
	var dx = animator.dx;
	var dy = animator.dy;
	var scale = animator.scale;

	g.fillRect(dx, dy, w * scale,h * scale);
	var proj = new Vector(0, 0);
	var val = new Vector(0, 0);
	g.lineWidth = .75;
	for (var i = 0; i < this.particles.length; i++) {
		var p = this.particles[i];
		if (!this.field.inBounds(p.x, p.y)) {
			p.age = -2;
			continue;
		}
		this.projection.project(p.x, p.y, proj);
		proj.x = proj.x * scale + dx;
		proj.y = proj.y * scale + dy;
		if (proj.x < 0 || proj.y < 0 || proj.x > w || proj.y > h) {
			p.age = -2;
		}
		if (p.oldX != -1) {
			var wind = this.field.getValue(p.x, p.y, val);
			var s = wind.length() / this.maxLength;
			var c = 90 + Math.round(350 * s); // was 400
			if (c > 255) {
				c = 255;
			} 
			g.strokeStyle = this.colors[c];
			g.beginPath();
			g.moveTo(proj.x, proj.y);
			g.lineTo(p.oldX, p.oldY);
			g.stroke();
	  }
		p.oldX = proj.x;
		p.oldY = proj.y;
	}
};

// please don't hate on this code too much.
// it's late and i'm tired.

var MotionDetails = function(div, callout, field, projection, animator) {
	$(callout).fadeOut();
	var moveTime = +new Date();
	var calloutOK = false;
	var currentlyShowing = false;
	var calloutX = 0;
	var calloutY = 0;
	var calloutHTML = '';
	var lastX = 0;
	var lastY = 0;

	function format(x) {
		x = Math.round(x * 10) / 10;
		var a1 = ~~x;
		var a2 = (~~(x * 10)) % 10;
		return a1 + '.' + a2;	
  } 

  function minutes(x) {	
		x = Math.round(x * 60) / 60;
		var degrees = ~~x;
		var m = ~~((x - degrees) * 60);
		return degrees + '&deg;&nbsp;' + (m == 0 ? '00' : m < 10 ? '0' + m : '' + m) + "'";
	}
	
	$(div).mouseleave(function() {
		moveTime = +new Date();
		calloutOK = false;
	});
	
	var pos = $(div).position();

	$(div).mousemove(function(e) {
		
		// TODO: REMOVE MAGIC CONSTANTS
		var x = e.pageX - this.offsetLeft - 60;
	  var y = e.pageY - this.offsetTop - 10;
		if (x == lastX && y == lastY) {
			return;
		}
		lastX = x;
		lastY = y;
		moveTime = +new Date();
		var scale = animator.scale;
		var dx = animator.dx;
		var dy = animator.dy;
		var mx = (x - dx) / scale;
		var my = (y - dy) / scale;
		var location = projection.invert(mx, my);
		var lat = location.y;
		var lon = location.x;
		var speed = 0;
		if (field.inBounds(lon, lat)) {
		  speed = field.getValue(lon, lat).length() / 1.15;
	  }
		calloutOK = !!speed;
		calloutHTML = '<div style="padding-bottom:5px"><b>' +
		              format(speed / 10)  + ' m/s</b> current speed<br></div>' +
		              minutes(lat)  + 'N, ' +
		              minutes(lon) + 'E<br>' +
									'click to zoom <br>';
	                  
		calloutY = (pos.top + y) + 'px';
		calloutX = (pos.left + x + 20) + 'px';
	});
	
	setInterval(function() {
		var timeSinceMove = +new Date() - moveTime;
		if (timeSinceMove > 200 && calloutOK) {
			if (!currentlyShowing) {
	  	  callout.innerHTML = calloutHTML;
				callout.style.left = calloutX;
				callout.style.top = calloutY;
				callout.style.visibility = 'visible';
				$(callout).fadeTo(400, 1);
				currentlyShowing = true;
			}
		} else if (currentlyShowing) {
	  	$(callout).fadeOut('fast');
			currentlyShowing = false;
		}
	}, 50);
};

/**
 * The cities array contains objects with properties: city, state, lat, lon, pop.
 *
 * @param {Array.<Object>} cities
 * @param {Object} canvas
 * @param {Object} projection
 */
var CityDisplay = function(cities, canvas, projection) {
	this.cities = cities;
	this.canvas = canvas;
	this.projection = projection;
	this.maxInView = 10;
	this.pad = 3;
	cities.sort(function(a, b) {
		return b.pop - a.pop;
	});
	for (var i = 0; i < this.cities.length; i++) {
		this.cities[i].alpha = 0;
	}
};

CityDisplay.prototype.endMove = function(animator) {
	for (var i = 0; i < this.cities.length; i++) {
		this.cities[i].alpha = 0;
	}
	this.move(animator);
}

CityDisplay.prototype.markCities = function(scale, dx, dy, alpha) {
	var spaceTaken = [];
	function collide(r1, r2) {
		return !(r1.x + r1.w < r2.x || r1.x > r2.x + r2.w ||
		         r1.y + r1.h < r2.y || r1.y > r2.y + r2.h);
	}
	function isFree(r) {
		for (var i = 0; i < spaceTaken.length; i++) {
			if (collide(r, spaceTaken[i])) {
				return false;
			}
		}
		return true;
	}
	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	var numInView = 0;
	var pad = this.pad;
	for (var i = 0; i < this.cities.length; i++) {
		var city = this.cities[i];
		var r = .075 * Math.pow(city.pop, .3);
		var v = this.projection.project(city.lon, city.lat);
		var x = v.x * scale + dx;
		var y = v.y * scale + dy;
		if (x < 0 || x > w || y < 0 || y > h) {
			continue;
		}
		var tx = x;
		var ty = y + 15;
	  var textSize = g.measureText(city.city);
		// check for collisions with previously drawn stuff.
		var dotArea = {
		  'x': x - r - pad,
		  'y': y - r - pad,
		  'w': 2 * (r + pad),
		  'h': 2 * (r + pad)
	  };
		var textArea = {
			'x': tx - textSize.width / 2 - pad,
			'y': ty - 15 - pad,
			'w': textSize.width + 2 * pad,
			'h': 15 + 2 * pad
		};
		if (!isFree(dotArea) || !isFree(textArea)) {
			continue;
		}
		spaceTaken.push(textArea);
		spaceTaken.push(dotArea);
		city.alpha += alpha;
		if (++numInView > this.maxInView) {
			break;
		}
	}

};


CityDisplay.prototype.move = function(animator) {
	for (var i = 0; i < this.cities.length; i++) {
	  this.cities[i].alpha = 0;
	}
	var dx = 0;
	var dy = 0
	var scale = 1;
	if (animator) {
		dx = animator.dx;
		dy = animator.dy;
		scale = animator.scale;
		if (animator.state == 'zoom') {

			var u = animator.zoomProgress;
			this.markCities(animator.scaleStart, animator.dxStart, animator.dyStart, 1 - u);
			this.markCities(animator.scaleTarget, animator.dxTarget, animator.dyTarget, u);
	  } else {
			//this.markCities(animator.scaleStart, animator.dxStart, animator.dyStart, 1);
			this.markCities(scale, dx, dy, 1);
		}
	} else {
		this.markCities(1, 0, 0, 1);
	}
	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	g.clearRect(0, 0, w, h);
	var pad = this.pad;
	for (var i = 0; i < this.cities.length; i++) {
		
		var city = this.cities[i];
		var alpha = Math.min(1, city.alpha);
		if (!alpha) {
			continue;
		}
		function check(val, name) {
			if (!val) {
				window.console.log(name + ' = ' + val);
			}
		}
		var r = .075 * Math.pow(city.pop, .3);
		var v = this.projection.project(city.lon, city.lat);
		var x = v.x * scale + dx;
		var y = v.y * scale + dy;

		if (x < 0 || x > w || y < 0 || y > h) {
			continue;
		}
		var tx = x;
		var ty = y + 15;

		g.beginPath();
    g.arc(x, y, r, 0, Math.PI*2, true); 
		g.closePath();
		g.fillStyle = 'rgba(255,255,255,' + alpha + ')';
		g.fill();
		g.strokeStyle = 'rgba(0,0,0,' + alpha + ')';
		g.stroke();
		g.fillStyle = 'rgba(255,255,255,' + .25 * alpha + ')';
		g.textAlign = 'center';
		g.font = '12px Verdana';

		for (var a = -2; a <= 2; a++) {
			for (var b = -2; b <= 2; b++) {
				g.fillText(city.city, tx + a, ty + b);
			}
		}
		g.fillStyle = 'rgba(0,0,0,' + alpha + ')';
		g.fillText(city.city, tx, ty);
	}
};
var cities = [
// I don't know why you're reading this boring file,
// but please note that "pop" isn't always a population number.
// The reasons are complicated, drop me a line if you really
// care a lot about it :-)
{
 city: 'Idaho Falls',
 state: 'Idaho',
 lat: 43.488673,
 lon: -112.03638,
 pop: 100000
},
{
 city: 'Amherst',
 state: 'Massachusetts',
 lat: 42.375623,
 lon: -72.51883,
 pop: 400000
},
{
 city: 'Lawrence',
 state: 'Kansas',
 lat: 38.96086,
 lon: -95.26454,
 pop: 100000
},
{
 city: 'Winchester',
 state: 'Massachusetts',
 lat: 42.45165,
 lon: -71.14659,
 pop: 100000
},
{
 city: 'Lexington',
 state: 'Massachusetts',
 lat: 42.44553,
 lon: -71.23065,
 pop: 100000
},
{
 city: 'Alliance',
 state: 'Nebraska',
 lat: 42.1001132,
 lon: -102.87543,
 pop: 100000
},
{
 city: 'Casper',
 state: 'Wyoming',
 lat: 42.83611,
 lon: -106.3452,
 pop: 100000
},
{
 city: 'Bangor',
 state: 'Maine',
 lat: 44.832518,
 lon: -68.790409,
 pop: 100000
},
{
 city: 'Grand Junction',
 state: 'Colorado',
 lat: 39.086035,
 lon: -108.56695,
 pop: 100000
},
{
 city: 'Rapid City',
 state: 'South Dakota',
 lat: 44.0700794,
 lon: -103.2178285,
 pop: 100000
},
{
 city: 'Bismarck',
 state: 'North Dakota',
 lat: 46.811371,
 lon: -100.768739,
 pop: 100000
},
{
 city: 'New York',
 state: 'New York',
 lat: 40.7143528,
 lon: -74.0059731,
 pop: 8175133
},
{
 city: 'Los Angeles',
 state: 'California',
 lat: 34.0522342,
 lon: -118.2436,
 pop: 3792621
},
{
 city: 'Chicago',
 state: 'Illinois',
 lat: 41.8781136,
 lon: -87.6297981,
 pop: 2695598
},
{
 city: 'Houston',
 state: 'Texas',
 lat: 29.7601927,
 lon: -95.369389,
 pop: 2099451
},
{
 city: 'Philadelphia',
 state: 'Pennsylvania',
 lat: 39.952335,
 lon: -75.163789,
 pop: 1526006
},
{
 city: 'Phoenix',
 state: 'Arizona',
 lat: 33.4483771,
 lon: -112.0740372,
 pop: 1445632
},
{
 city: 'San Antonio',
 state: 'Texas',
 lat: 29.4241219,
 lon: -98.4936281,
 pop: 1327407
},
{
 city: 'San Diego',
 state: 'California',
 lat: 32.7153292,
 lon: -117.157255,
 pop: 1307402
},
{
 city: 'Dallas',
 state: 'Texas',
 lat: 32.802955,
 lon: -96.769923,
 pop: 1197816
},
{
 city: 'San Jose',
 state: 'California',
 lat: 37.3393857,
 lon: -121.8949554,
 pop: 945942
},
{
 city: 'Jacksonville',
 state: 'Florida',
 lat: 30.3321838,
 lon: -81.655650,
 pop: 821784
},
{
 city: 'Indianapolis',
 state: 'Indiana',
 lat: 39.7685155,
 lon: -86.1580735,
 pop: 820445
},
{
 city: 'San Francisco',
 state: 'California',
 lat: 37.7749295,
 lon: -122.4194155,
 pop: 805235
},
{
 city: 'Austin',
 state: 'Texas',
 lat: 30.267153,
 lon: -97.7430607,
 pop: 790390
},
{
 city: 'Columbus',
 state: 'Ohio',
 lat: 39.9611755,
 lon: -82.9987942,
 pop: 1787033.1
},
{
 city: 'Fort Worth',
 state: 'Texas',
 lat: 32.725409,
 lon: -97.3208495,
 pop: 741206
},
{
 city: 'Charlotte',
 state: 'North Carolina',
 lat: 35.2270869,
 lon: -80.8431266,
 pop: 731424
},
{
 city: 'Detroit',
 state: 'Michigan',
 lat: 42.331427,
 lon: -83.0457538,
 pop: 713777
},
{
 city: 'El Paso',
 state: 'Texas',
 lat: 31.7587198,
 lon: -106.4869314,
 pop: 649121
},
{
 city: 'Memphis',
 state: 'Tennessee',
 lat: 35.1495343,
 lon: -90.0489801,
 pop: 646889
},
{
 city: 'Baltimore',
 state: 'Maryland',
 lat: 39.2903848,
 lon: -76.6121893,
 pop: 620961
},
{
 city: 'Boston',
 state: 'Massachusetts',
 lat: 42.3584308,
 lon: -71.0597732,
 pop: 617594
},
{
 city: 'Seattle',
 state: 'Washington',
 lat: 47.6062095,
 lon: -122.3320708,
 pop: 1608660.1
},
{
 city: 'Washington',
 state: 'District of Columbia',
 lat: 38.8951118,
 lon: -77.0363658,
 pop: 601723
},
{
 city: 'Nashville',
 state: 'Tennessee',
 lat: 36.1658899,
 lon: -86.7844432,
 pop: 601222
},
{
 city: 'Denver',
 state: 'Colorado',
 lat: 39.7391536,
 lon: -104.9847034,
 pop: 1200000
},
{
 city: 'Louisville',
 state: 'Kentucky',
 lat: 38.2526647,
 lon: -85.75845571,
 pop: 597337
},
{
 city: 'Milwaukee',
 state: 'Wisconsin',
 lat: 43.0389025,
 lon: -87.90647363,
 pop: 594833
},
{
 city: 'Portland',
 state: 'Oregon',
 lat: 45.5234515,
 lon: -122.6762071,
 pop: 583776
},
{
 city: 'Las Vegas',
 state: 'Nevada',
 lat: 36.114646,
 lon: -115.17281601,
 pop: 583756
},
{
 city: 'Oklahoma City',
 state: 'Oklahoma',
 lat: 35.4675602,
 lon: -97.51642759,
 pop: 579999
},
{
 city: 'Albuquerque',
 state: 'New Mexico',
 lat: 35.0844909,
 lon: -106.6511367,
 pop: 545852
},
{
 city: 'Tucson',
 state: 'Arizona',
 lat: 32.2217429,
 lon: -110.92647897,
 pop: 520116
},
{
 city: 'Fresno',
 state: 'California',
 lat: 36.7477272,
 lon: -119.7723661,
 pop: 494665
},
{
 city: 'Sacramento',
 state: 'California',
 lat: 38.5815719,
 lon: -121.49439961,
 pop: 466488
},
{
 city: 'Long Beach',
 state: 'California',
 lat: 33.8041667,
 lon: -118.15805561,
 pop: 462257
},
{
 city: 'Kansas City',
 state: 'Missouri',
 lat: 39.0997265,
 lon: -94.57856671,
 pop: 459787
},
{
 city: 'Mesa',
 state: 'Arizona',
 lat: 33.4151843,
 lon: -111.8314724,
 pop: 439041
},
{
 city: 'Virginia Beach',
 state: 'Virginia',
 lat: 36.8529263,
 lon: -75.97798499,
 pop: 437994
},
{
 city: 'Atlanta',
 state: 'Georgia',
 lat: 33.7489954,
 lon: -84.3879824,
 pop: 420003
},
{
 city: 'Colorado Springs',
 state: 'Colorado',
 lat: 38.8338816,
 lon: -104.8213634,
 pop: 416427
},
{
 city: 'Omaha',
 state: 'Nebraska',
 lat: 41.2523634,
 lon: -95.99798827,
 pop: 408958
},
{
 city: 'Raleigh',
 state: 'North Carolina',
 lat: 35.772096,
 lon: -78.63861452,
 pop: 403892
},
{
 city: 'Miami',
 state: 'Florida',
 lat: 25.7889689,
 lon: -80.22643928,
 pop: 399457
},
{
 city: 'Cleveland',
 state: 'Ohio',
 lat: 41.4994954,
 lon: -81.6954088,
 pop: 396815
},
{
 city: 'Tulsa',
 state: 'Oklahoma',
 lat: 36.1539816,
 lon: -95.992775,
 pop: 391906
},
{
 city: 'Oakland',
 state: 'California',
 lat: 37.8043637,
 lon: -122.2711137,
 pop: 390724
},
{
 city: 'Minneapolis',
 state: 'Minnesota',
 lat: 44.983334,
 lon: -93.26666998,
 pop: 382578
},
{
 city: 'Wichita',
 state: 'Kansas',
 lat: 37.6922222,
 lon: -97.33722219,
 pop: 382368
},
{
 city: 'Arlington',
 state: 'Texas',
 lat: 32.735687,
 lon: -97.10806557,
 pop: 365438
},
{
 city: 'Bakersfield',
 state: 'California',
 lat: 35.3732921,
 lon: -119.01871249,
 pop: 347483
},
{
 city: 'New Orleans',
 state: 'Louisiana',
 lat: 29.95106579,
 lon: -90.0715323,
 pop: 343829
},
{
 city: 'Honolulu',
 state: 'Hawaii',
 lat: 21.3069444,
 lon: -157.85833333,
 pop: 337256
},
{
 city: 'Anaheim',
 state: 'California',
 lat: 33.8352932,
 lon: -117.91450359,
 pop: 336265
},
{
 city: 'Tampa',
 state: 'Florida',
 lat: 27.950575,
 lon: -82.45717762,
 pop: 335709
},
{
 city: 'Aurora',
 state: 'Colorado',
 lat: 39.7294319,
 lon: -104.83191953,
 pop: 325078
},
{
 city: 'Santa Ana',
 state: 'California',
 lat: 33.7455731,
 lon: -117.86783377,
 pop: 324528
},
{
 city: 'Saint Louis',
 state: 'Missouri',
 lat: 38.6270025,
 lon: -90.1994042,
 pop: 319294
},
{
 city: 'Pittsburgh',
 state: 'Pennsylvania',
 lat: 40.44062479,
 lon: -79.99588642,
 pop: 305704
},
{
 city: 'Corpus Christi',
 state: 'Texas',
 lat: 27.8005828,
 lon: -97.39638102,
 pop: 305215
},
{
 city: 'Riverside',
 state: 'California',
 lat: 33.9533487,
 lon: -117.3961564,
 pop: 303871
},
{
 city: 'Cincinnati',
 state: 'Ohio',
 lat: 39.1031182,
 lon: -84.51201963,
 pop: 296943
},
{
 city: 'Lexington',
 state: 'Kentucky',
 lat: 38.0405837,
 lon: -84.50371643,
 pop: 295803
},
{
 city: 'Anchorage',
 state: 'Alaska',
 lat: 61.2180556,
 lon: -149.90027783,
 pop: 291826
},
{
 city: 'Stockton',
 state: 'California',
 lat: 37.9577016,
 lon: -121.29077961,
 pop: 291707
},
{
 city: 'Toledo',
 state: 'Ohio',
 lat: 41.6639383,
 lon: -83.55521198,
 pop: 287208
},
{
 city: 'Saint Paul',
 state: 'Minnesota',
 lat: 44.95416669,
 lon: -93.1138889,
 pop: 285068
},
{
 city: 'Newark',
 state: 'New Jersey',
 lat: 40.735657,
 lon: -74.1723667,
 pop: 277140
},
{
 city: 'Greensboro',
 state: 'North Carolina',
 lat: 36.0726354,
 lon: -79.79197541,
 pop: 269666
},
{
 city: 'Buffalo',
 state: 'New York',
 lat: 42.88644679,
 lon: -78.8783689,
 pop: 261310
},
{
 city: 'Plano',
 state: 'Texas',
 lat: 33.0198431,
 lon: -96.69888558,
 pop: 259841
},
{
 city: 'Lincoln',
 state: 'Nebraska',
 lat: 40.806862,
 lon: -96.68167903,
 pop: 258379
},
{
 city: 'Henderson',
 state: 'Nevada',
 lat: 36.0395247,
 lon: -114.9817213,
 pop: 257729
},
{
 city: 'Fort Wayne',
 state: 'Indiana',
 lat: 41.079273,
 lon: -85.13935129,
 pop: 253691
},
{
 city: 'Jersey City',
 state: 'New Jersey',
 lat: 40.72815749,
 lon: -74.07764172,
 pop: 247597
},
{
 city: 'Saint Petersburg',
 state: 'Florida',
 lat: 27.7730556,
 lon: -82.63,
 pop: 244769
},
{
 city: 'Chula Vista',
 state: 'California',
 lat: 32.6400541,
 lon: -117.08419552,
 pop: 243916
},
{
 city: 'Norfolk',
 state: 'Virginia',
 lat: 36.8507689,
 lon: -76.2858726,
 pop: 242803
},
{
 city: 'Orlando',
 state: 'Florida',
 lat: 28.5383355,
 lon: -81.37923649,
 pop: 238300
},
{
 city: 'Chandler',
 state: 'Arizona',
 lat: 33.3061605,
 lon: -111.84125019,
 pop: 236123
},
{
 city: 'Laredo',
 state: 'Texas',
 lat: 27.506407,
 lon: -99.50754212,
 pop: 236091
},
{
 city: 'Madison',
 state: 'Wisconsin',
 lat: 43.0730517,
 lon: -89.40123019,
 pop: 233209
},
{
 city: 'Winston-Salem',
 state: 'North Carolina',
 lat: 36.09985959,
 lon: -80.244216,
 pop: 229617
},
{
 city: 'Lubbock',
 state: 'Texas',
 lat: 33.5778631,
 lon: -101.8551665,
 pop: 229573
},
{
 city: 'Baton Rouge',
 state: 'Louisiana',
 lat: 30.4582829,
 lon: -91.1403196,
 pop: 229493
},
{
 city: 'Durham',
 state: 'North Carolina',
 lat: 35.9940329,
 lon: -78.898619,
 pop: 228330
},
{
 city: 'Garland',
 state: 'Texas',
 lat: 32.912624,
 lon: -96.63888327,
 pop: 226876
},
{
 city: 'Glendale',
 state: 'Arizona',
 lat: 33.5386523,
 lon: -112.18598658,
 pop: 226721
},
{
 city: 'Reno',
 state: 'Nevada',
 lat: 39.5296329,
 lon: -119.8138027,
 pop: 225221
},
{
 city: 'Hialeah',
 state: 'Florida',
 lat: 25.8575963,
 lon: -80.27810567,
 pop: 224669
},
{
 city: 'Chesapeake [e]',
 state: 'Virginia',
 lat: 36.7682088,
 lon: -76.28749273,
 pop: 222209
},
{
 city: 'Scottsdale',
 state: 'Arizona',
 lat: 33.4941704,
 lon: -111.9260519,
 pop: 217385
},
{
 city: 'North Las Vegas',
 state: 'Nevada',
 lat: 36.1988592,
 lon: -115.11750131,
 pop: 216961
},
{
 city: 'Irving',
 state: 'Texas',
 lat: 32.8140177,
 lon: -96.9488945,
 pop: 216290
},
{
 city: 'Fremont',
 state: 'California',
 lat: 37.5482697,
 lon: -121.98857191,
 pop: 214089
},
{
 city: 'Irvine',
 state: 'California',
 lat: 33.6839473,
 lon: -117.79469418,
 pop: 212375
},
{
 city: 'Birmingham',
 state: 'Alabama',
 lat: 33.5206608,
 lon: -86.80248998,
 pop: 212237
},
{
 city: 'Rochester',
 state: 'New York',
 lat: 43.16103,
 lon: -77.6109219,
 pop: 210565
},
{
 city: 'San Bernardino',
 state: 'California',
 lat: 34.1083449,
 lon: -117.28976523,
 pop: 209924
},
{
 city: 'Spokane',
 state: 'Washington',
 lat: 47.6587802,
 lon: -117.4260466,
 pop: 208916
},
{
 city: 'Gilbert',
 state: 'Arizona',
 lat: 33.3528264,
 lon: -111.78902703,
 pop: 208453
},
{
 city: 'Arlington',
 state: 'Virginia',
 lat: 38.8799697,
 lon: -77.1067698,
 pop: 207627
},
{
 city: 'Montgomery',
 state: 'Alabama',
 lat: 32.3668052,
 lon: -86.29996891,
 pop: 205764
},
{
 city: 'Boise',
 state: 'Idaho',
 lat: 43.612631,
 lon: -116.21107599,
 pop: 205671
},
{
 city: 'Richmond',
 state: 'Virginia',
 lat: 37.5407246,
 lon: -77.4360481,
 pop: 204214
},
{
 city: 'Des Moines',
 state: 'Iowa',
 lat: 41.6005448,
 lon: -93.60910637,
 pop: 203433
},
{
 city: 'Modesto',
 state: 'California',
 lat: 37.63909719,
 lon: -120.99687817,
 pop: 201165
},
{
 city: 'Fayetteville',
 state: 'North Carolina',
 lat: 35.0526641,
 lon: -78.87835849,
 pop: 200654
},
{
 city: 'Shreveport',
 state: 'Louisiana',
 lat: 32.5251516,
 lon: -93.75017888,
 pop: 199311
},
{
 city: 'Akron',
 state: 'Ohio',
 lat: 41.0814447,
 lon: -81.5190053,
 pop: 199110
},
{
 city: 'Tacoma',
 state: 'Washington',
 lat: 47.2528768,
 lon: -122.44429059,
 pop: 198397
},
{
 city: 'Aurora',
 state: 'Illinois',
 lat: 41.7605849,
 lon: -88.32007154,
 pop: 197899
},
{
 city: 'Oxnard',
 state: 'California',
 lat: 34.1975048,
 lon: -119.17705163,
 pop: 197899
},
{
 city: 'Fontana',
 state: 'California',
 lat: 34.0922335,
 lon: -117.435048,
 pop: 196069
},
{
 city: 'Yonkers',
 state: 'New York',
 lat: 40.9312099,
 lon: -73.89874689,
 pop: 195976
},
{
 city: 'Augusta',
 state: 'Georgia',
 lat: 33.474246,
 lon: -82.00967003,
 pop: 195844
},
{
 city: 'Mobile',
 state: 'Alabama',
 lat: 30.6943566,
 lon: -88.0430541,
 pop: 195111
},
{
 city: 'Little Rock',
 state: 'Arkansas',
 lat: 34.7464809,
 lon: -92.28959477,
 pop: 193524
},
{
 city: 'Moreno Valley',
 state: 'California',
 lat: 33.9424658,
 lon: -117.22967168,
 pop: 193365
},
{
 city: 'Glendale',
 state: 'California',
 lat: 34.1425078,
 lon: -118.25507503,
 pop: 191719
},
{
 city: 'Amarillo',
 state: 'Texas',
 lat: 35.2219971,
 lon: -101.83129688,
 pop: 190695
},
{
 city: 'Huntington Beach',
 state: 'California',
 lat: 33.660297,
 lon: -117.99922652,
 pop: 189992
},
{
 city: 'Columbus',
 state: 'Georgia',
 lat: 32.4609764,
 lon: -84.98770937,
 pop: 189885
},
{
 city: 'Grand Rapids',
 state: 'Michigan',
 lat: 42.9633599,
 lon: -85.66808633,
 pop: 188040
},
{
 city: 'Salt Lake City',
 state: 'Utah',
 lat: 40.7607793,
 lon: -111.89104739,
 pop: 186440
},
{
 city: 'Tallahassee',
 state: 'Florida',
 lat: 30.4382559,
 lon: -84.28073288,
 pop: 181376
},
{
 city: 'Worcester',
 state: 'Massachusetts',
 lat: 42.2625932,
 lon: -71.8022934,
 pop: 181045
},
{
 city: 'Newport News',
 state: 'Virginia',
 lat: 37.0870821,
 lon: -76.47301217,
 pop: 180719
},
{
 city: 'Huntsville',
 state: 'Alabama',
 lat: 34.7303688,
 lon: -86.58610367,
 pop: 180105
},
{
 city: 'Knoxville',
 state: 'Tennessee',
 lat: 35.9606384,
 lon: -83.92073921,
 pop: 178874
},
{
 city: 'Providence',
 state: 'Rhode Island',
 lat: 41.8239891,
 lon: -71.41283429,
 pop: 178042
},
{
 city: 'Santa Clarita',
 state: 'California',
 lat: 34.3916641,
 lon: -118.54258603,
 pop: 176320
},
{
 city: 'Grand Prairie',
 state: 'Texas',
 lat: 32.7459645,
 lon: -96.99778459,
 pop: 175396
},
{
 city: 'Brownsville',
 state: 'Texas',
 lat: 25.9017472,
 lon: -97.4974838,
 pop: 175023
},
{
 city: 'Jackson',
 state: 'Mississippi',
 lat: 32.2987573,
 lon: -90.18481028,
 pop: 173514
},
{
 city: 'Overland Park',
 state: 'Kansas',
 lat: 38.9822282,
 lon: -94.6707917,
 pop: 173372
},
{
 city: 'Garden Grove',
 state: 'California',
 lat: 33.7739053,
 lon: -117.94144773,
 pop: 170883
},
{
 city: 'Santa Rosa',
 state: 'California',
 lat: 38.4404674,
 lon: -122.71443137,
 pop: 167815
},
{
 city: 'Chattanooga',
 state: 'Tennessee',
 lat: 35.0456297,
 lon: -85.30968008,
 pop: 167674
},
{
 city: 'Oceanside',
 state: 'California',
 lat: 33.1958696,
 lon: -117.37948343,
 pop: 167086
},
{
 city: 'Fort Lauderdale',
 state: 'Florida',
 lat: 26.1223084,
 lon: -80.1433786,
 pop: 165521
},
{
 city: 'Rancho Cucamonga',
 state: 'California',
 lat: 34.10639889,
 lon: -117.5931084,
 pop: 165269
},
{
 city: 'Port Saint Lucie',
 state: 'Florida',
 lat: 27.2758333,
 lon: -80.35500002,
 pop: 164603
},
{
 city: 'Ontario',
 state: 'California',
 lat: 34.0633443,
 lon: -117.65088763,
 pop: 163924
},
{
 city: 'Vancouver',
 state: 'Washington',
 lat: 45.6387281,
 lon: -122.66148609,
 pop: 161791
},
{
 city: 'Tempe',
 state: 'Arizona',
 lat: 33.4255104,
 lon: -111.94000542,
 pop: 161719
},
{
 city: 'Springfield',
 state: 'Missouri',
 lat: 37.2089572,
 lon: -93.29229889,
 pop: 159498
},
{
 city: 'Lancaster',
 state: 'California',
 lat: 34.6867846,
 lon: -118.15416317,
 pop: 156633
},
{
 city: 'Eugene',
 state: 'Oregon',
 lat: 44.0520691,
 lon: -123.08675361,
 pop: 156185
},
{
 city: 'Pembroke Pines',
 state: 'Florida',
 lat: 26.0122378,
 lon: -80.31522331,
 pop: 154750
},
{
 city: 'Salem',
 state: 'Oregon',
 lat: 44.9428975,
 lon: -123.03509632,
 pop: 154637
},
{
 city: 'Cape Coral',
 state: 'Florida',
 lat: 26.5628537,
 lon: -81.9495331,
 pop: 154305
},
{
 city: 'Peoria',
 state: 'Arizona',
 lat: 33.5805955,
 lon: -112.23737791,
 pop: 154065
},
{
 city: 'Sioux Falls',
 state: 'South Dakota',
 lat: 43.5499749,
 lon: -96.70032702,
 pop: 153888
},
{
 city: 'Springfield',
 state: 'Massachusetts',
 lat: 42.1014831,
 lon: -72.589811,
 pop: 153060
},
{
 city: 'Elk Grove',
 state: 'California',
 lat: 38.4087993,
 lon: -121.37161777,
 pop: 153015
},
{
 city: 'Rockford',
 state: 'Illinois',
 lat: 42.2711311,
 lon: -89.0939952,
 pop: 152871
},
{
 city: 'Palmdale',
 state: 'California',
 lat: 34.5794343,
 lon: -118.11646127,
 pop: 152750
},
{
 city: 'Corona',
 state: 'California',
 lat: 33.8752935,
 lon: -117.56643838,
 pop: 152374
},
{
 city: 'Salinas',
 state: 'California',
 lat: 36.6777372,
 lon: -121.65550127,
 pop: 150441
},
{
 city: 'Pomona',
 state: 'California',
 lat: 34.0552267,
 lon: -117.75230479,
 pop: 149058
},
{
 city: 'Pasadena',
 state: 'Texas',
 lat: 29.6910625,
 lon: -95.2091006,
 pop: 149043
},
{
 city: 'Joliet',
 state: 'Illinois',
 lat: 41.525031,
 lon: -88.08172507,
 pop: 147433
},
{
 city: 'Paterson',
 state: 'New Jersey',
 lat: 40.9167654,
 lon: -74.17181099,
 pop: 146199
},
{
 city: 'Kansas City',
 state: 'Kansas',
 lat: 39.114053,
 lon: -94.6274636,
 pop: 145786
},
{
 city: 'Torrance',
 state: 'California',
 lat: 33.8358492,
 lon: -118.34062879,
 pop: 145438
},
{
 city: 'Syracuse',
 state: 'New York',
 lat: 43.0481221,
 lon: -76.14742438,
 pop: 145170
},
{
 city: 'Bridgeport',
 state: 'Connecticut',
 lat: 41.1865478,
 lon: -73.19517669,
 pop: 144229
},
{
 city: 'Hayward',
 state: 'California',
 lat: 37.6688205,
 lon: -122.0807964,
 pop: 144186
},
{
 city: 'Fort Collins',
 state: 'Colorado',
 lat: 40.5852602,
 lon: -105.08442302,
 pop: 143986
},
{
 city: 'Escondido',
 state: 'California',
 lat: 33.1192068,
 lon: -117.08642097,
 pop: 143911
},
{
 city: 'Lakewood',
 state: 'Colorado',
 lat: 39.7047095,
 lon: -105.08137342,
 pop: 142980
},
{
 city: 'Naperville',
 state: 'Illinois',
 lat: 41.7858629,
 lon: -88.14728931,
 pop: 141853
},
{
 city: 'Dayton',
 state: 'Ohio',
 lat: 39.7589478,
 lon: -84.19160691,
 pop: 141527
},
{
 city: 'Hollywood',
 state: 'Florida',
 lat: 26.0112014,
 lon: -80.14949008,
 pop: 140768
},
{
 city: 'Sunnyvale',
 state: 'California',
 lat: 37.36883,
 lon: -122.0363496,
 pop: 140081
},
{
 city: 'Alexandria',
 state: 'Virginia',
 lat: 38.8048355,
 lon: -77.04692137,
 pop: 139966
},
{
 city: 'Mesquite',
 state: 'Texas',
 lat: 32.76679551,
 lon: -96.5991593,
 pop: 139824
},
{
 city: 'Hampton',
 state: 'Virginia',
 lat: 37.0298687,
 lon: -76.34522179,
 pop: 137436
},
{
 city: 'Pasadena',
 state: 'California',
 lat: 34.1477849,
 lon: -118.14451551,
 pop: 137122
},
{
 city: 'Orange',
 state: 'California',
 lat: 33.7877944,
 lon: -117.85311189,
 pop: 136416
},
{
 city: 'Savannah',
 state: 'Georgia',
 lat: 32.0835407,
 lon: -81.09983418,
 pop: 136286
},
{
 city: 'Cary',
 state: 'North Carolina',
 lat: 35.79154,
 lon: -78.78111693,
 pop: 135234
},
{
 city: 'Fullerton',
 state: 'California',
 lat: 33.8702923,
 lon: -117.92533801,
 pop: 135161
},
{
 city: 'Warren',
 state: 'Michigan',
 lat: 42.49299999,
 lon: -83.02819698,
 pop: 134056
},
{
 city: 'Clarksville',
 state: 'Tennessee',
 lat: 36.5297706,
 lon: -87.35945279,
 pop: 132929
},
{
 city: 'McKinney',
 state: 'Texas',
 lat: 33.1972465,
 lon: -96.63978221,
 pop: 131117
},
{
 city: 'McAllen',
 state: 'Texas',
 lat: 26.2034071,
 lon: -98.23001236,
 pop: 129877
},
{
 city: 'New Haven',
 state: 'Connecticut',
 lat: 41.3081527,
 lon: -72.92815769,
 pop: 129779
},
{
 city: 'Sterling Heights',
 state: 'Michigan',
 lat: 42.5803122,
 lon: -83.03020328,
 pop: 129699
},
{
 city: 'West Valley City',
 state: 'Utah',
 lat: 40.6916132,
 lon: -112.00105009,
 pop: 129480
},
{
 city: 'Columbia',
 state: 'South Carolina',
 lat: 34.0007104,
 lon: -81.03481442,
 pop: 129272
},
{
 city: 'Killeen',
 state: 'Texas',
 lat: 31.1171194,
 lon: -97.72779589,
 pop: 127921
},
{
 city: 'Topeka',
 state: 'Kansas',
 lat: 39.0558235,
 lon: -95.68901847,
 pop: 127473
},
{
 city: 'Thousand Oaks',
 state: 'California',
 lat: 34.1705609,
 lon: -118.83759371,
 pop: 126683
},
{
 city: 'Cedar Rapids',
 state: 'Iowa',
 lat: 41.9778795,
 lon: -91.66562323,
 pop: 126326
},
{
 city: 'Olathe',
 state: 'Kansas',
 lat: 38.8813958,
 lon: -94.81912848,
 pop: 125872
},
{
 city: 'Elizabeth',
 state: 'New Jersey',
 lat: 40.6639916,
 lon: -74.2107006,
 pop: 124969
},
{
 city: 'Waco',
 state: 'Texas',
 lat: 31.549333,
 lon: -97.14666953,
 pop: 124805
},
{
 city: 'Hartford',
 state: 'Connecticut',
 lat: 41.76371109,
 lon: -72.68509318,
 pop: 124775
},
{
 city: 'Visalia',
 state: 'California',
 lat: 36.3302284,
 lon: -119.2920585,
 pop: 124442
},
{
 city: 'Gainesville',
 state: 'Florida',
 lat: 29.6516344,
 lon: -82.32482616,
 pop: 124354
},
{
 city: 'Simi Valley',
 state: 'California',
 lat: 34.2694474,
 lon: -118.78148198,
 pop: 124237
},
{
 city: 'Stamford',
 state: 'Connecticut',
 lat: 41.0534302,
 lon: -73.5387341,
 pop: 122643
},
{
 city: 'Bellevue',
 state: 'Washington',
 lat: 47.610377,
 lon: -122.2006786,
 pop: 122363
},
{
 city: 'Concord',
 state: 'California',
 lat: 37.9779776,
 lon: -122.0310733,
 pop: 122067
},
{
 city: 'Miramar',
 state: 'Florida',
 lat: 25.9756704,
 lon: -80.2867501,
 pop: 122041
},
{
 city: 'Coral Springs',
 state: 'Florida',
 lat: 26.271192,
 lon: -80.27060442,
 pop: 121096
},
{
 city: 'Lafayette',
 state: 'Louisiana',
 lat: 30.2240897,
 lon: -92.01984273,
 pop: 120623
},
{
 city: 'Charleston',
 state: 'South Carolina',
 lat: 32.7765656,
 lon: -79.93092158,
 pop: 120083
},
{
 city: 'Carrollton',
 state: 'Texas',
 lat: 32.9756415,
 lon: -96.88996359,
 pop: 119097
},
{
 city: 'Roseville',
 state: 'California',
 lat: 38.7521235,
 lon: -121.28800593,
 pop: 118788
},
{
 city: 'Thornton',
 state: 'Colorado',
 lat: 39.8680412,
 lon: -104.97192431,
 pop: 118772
},
{
 city: 'Beaumont',
 state: 'Texas',
 lat: 30.080174,
 lon: -94.12655618,
 pop: 118296
},
{
 city: 'Allentown',
 state: 'Pennsylvania',
 lat: 40.6084305,
 lon: -75.49018331,
 pop: 118032
},
{
 city: 'Surprise',
 state: 'Arizona',
 lat: 33.639099,
 lon: -112.39575762,
 pop: 117517
},
{
 city: 'Evansville',
 state: 'Indiana',
 lat: 37.9715592,
 lon: -87.57108978,
 pop: 117429
},
{
 city: 'Abilene',
 state: 'Texas',
 lat: 32.4487364,
 lon: -99.73314392,
 pop: 117063
},
{
 city: 'Frisco',
 state: 'Texas',
 lat: 33.1506744,
 lon: -96.82361159,
 pop: 116989
},
{
 city: 'Independence',
 state: 'Missouri',
 lat: 39.0911161,
 lon: -94.4155068,
 pop: 116830
},
{
 city: 'Santa Clara',
 state: 'California',
 lat: 37.3541079,
 lon: -121.95523558,
 pop: 116468
},
{
 city: 'Springfield',
 state: 'Illinois',
 lat: 39.78172131,
 lon: -89.65014812,
 pop: 116250
},
{
 city: 'Vallejo',
 state: 'California',
 lat: 38.1040864,
 lon: -122.2566367,
 pop: 115942
},
{
 city: 'Victorville',
 state: 'California',
 lat: 34.5361067,
 lon: -117.2911565,
 pop: 115903
},
{
 city: 'Athens',
 state: 'Georgia',
 lat: 33.955802,
 lon: -83.38236561,
 pop: 115452
},
{
 city: 'Peoria',
 state: 'Illinois',
 lat: 40.6936488,
 lon: -89.58898641,
 pop: 115007
},
{
 city: 'Lansing',
 state: 'Michigan',
 lat: 42.732535,
 lon: -84.55553471,
 pop: 114297
},
{
 city: 'Ann Arbor',
 state: 'Michigan',
 lat: 42.2808256,
 lon: -83.74303782,
 pop: 113934
},
{
 city: 'El Monte',
 state: 'California',
 lat: 34.0686206,
 lon: -118.02756667,
 pop: 113475
},
{
 city: 'Denton',
 state: 'Texas',
 lat: 33.2148412,
 lon: -97.13306829,
 pop: 113383
},
{
 city: 'Berkeley',
 state: 'California',
 lat: 37.8715926,
 lon: -122.27274698,
 pop: 112580
},
{
 city: 'Provo',
 state: 'Utah',
 lat: 40.2338438,
 lon: -111.65853372,
 pop: 112488
},
{
 city: 'Downey',
 state: 'California',
 lat: 33.94001431,
 lon: -118.1325688,
 pop: 111772
},
{
 city: 'Midland',
 state: 'Texas',
 lat: 31.9973456,
 lon: -102.07791459,
 pop: 111147
},
{
 city: 'Norman',
 state: 'Oklahoma',
 lat: 35.2225668,
 lon: -97.4394777,
 pop: 110925
},
{
 city: 'Waterbury',
 state: 'Connecticut',
 lat: 41.5581525,
 lon: -73.05149648,
 pop: 110366
},
{
 city: 'Costa Mesa',
 state: 'California',
 lat: 33.6411316,
 lon: -117.9186689,
 pop: 109960
},
{
 city: 'Inglewood',
 state: 'California',
 lat: 33.9616801,
 lon: -118.35313108,
 pop: 109673
},
{
 city: 'Manchester',
 state: 'New Hampshire',
 lat: 42.9956397,
 lon: -71.45478907,
 pop: 109565
},
{
 city: 'Murfreesboro',
 state: 'Tennessee',
 lat: 35.8456213,
 lon: -86.39026999,
 pop: 108755
},
{
 city: 'Columbia',
 state: 'Missouri',
 lat: 38.9517053,
 lon: -92.33407237,
 pop: 108500
},
{
 city: 'Elgin',
 state: 'Illinois',
 lat: 42.0372487,
 lon: -88.28118948,
 pop: 108188
},
{
 city: 'Clearwater',
 state: 'Florida',
 lat: 27.9658533,
 lon: -82.8001026,
 pop: 107685
},
{
 city: 'Miami Gardens',
 state: 'Florida',
 lat: 25.9420377,
 lon: -80.24560451,
 pop: 107167
},
{
 city: 'Rochester',
 state: 'Minnesota',
 lat: 44.0216306,
 lon: -92.46989919,
 pop: 106769
},
{
 city: 'Pueblo',
 state: 'Colorado',
 lat: 38.2544472,
 lon: -104.6091409,
 pop: 106595
},
{
 city: 'Lowell',
 state: 'Massachusetts',
 lat: 42.6334247,
 lon: -71.3161718,
 pop: 106519
},
{
 city: 'Wilmington',
 state: 'North Carolina',
 lat: 34.2257255,
 lon: -77.94471023,
 pop: 106476
},
{
 city: 'Arvada',
 state: 'Colorado',
 lat: 39.8027644,
 lon: -105.0874842,
 pop: 106433
},
{
 city: 'San Buenaventura (Ventura)',
 state: 'California',
 lat: 34.2746405,
 lon: -119.22900528,
 pop: 106433
},
{
 city: 'Westminster',
 state: 'Colorado',
 lat: 39.8366528,
 lon: -105.0372046,
 pop: 106114
},
{
 city: 'West Covina',
 state: 'California',
 lat: 34.0686208,
 lon: -117.9389526,
 pop: 106098
},
{
 city: 'Gresham',
 state: 'Oregon',
 lat: 45.5001357,
 lon: -122.43020132,
 pop: 105594
},
{
 city: 'Fargo',
 state: 'North Dakota',
 lat: 46.8771863,
 lon: -96.78980338,
 pop: 105549
},
{
 city: 'Norwalk',
 state: 'California',
 lat: 33.9022367,
 lon: -118.08173299,
 pop: 105549
},
{
 city: 'Carlsbad',
 state: 'California',
 lat: 33.1580933,
 lon: -117.35059394,
 pop: 105328
},
{
 city: 'Fairfield',
 state: 'California',
 lat: 38.24935809,
 lon: -122.0399663,
 pop: 105321
},
{
 city: 'Cambridge',
 state: 'Massachusetts',
 lat: 42.3736158,
 lon: -71.1097335,
 pop: 105162
},
{
 city: 'Wichita Falls',
 state: 'Texas',
 lat: 33.9137085,
 lon: -98.4933873,
 pop: 104553
},
{
 city: 'Billings',
 state: 'Montana',
 lat: 45.7832856,
 lon: -108.5006904,
 pop: 104170
},
{
 city: 'Green Bay',
 state: 'Wisconsin',
 lat: 44.51915899,
 lon: -88.01982597,
 pop: 104057
},
{
 city: 'West Jordan',
 state: 'Utah',
 lat: 40.6096698,
 lon: -111.93910311,
 pop: 103712
},
{
 city: 'Richmond',
 state: 'California',
 lat: 37.9357576,
 lon: -122.34774859,
 pop: 103701
},
{
 city: 'Murrieta',
 state: 'California',
 lat: 33.5539143,
 lon: -117.21392321,
 pop: 103466
},
{
 city: 'Burbank',
 state: 'California',
 lat: 34.1808392,
 lon: -118.30896612,
 pop: 103340
},
{
 city: 'Palm Bay',
 state: 'Florida',
 lat: 28.0344621,
 lon: -80.58866462,
 pop: 103190
},
{
 city: 'Everett',
 state: 'Washington',
 lat: 47.9789848,
 lon: -122.2020794,
 pop: 103019
},
{
 city: 'Flint',
 state: 'Michigan',
 lat: 43.0125274,
 lon: -83.68745619,
 pop: 102434
},
{
 city: 'Antioch',
 state: 'California',
 lat: 38.0049214,
 lon: -121.805789,
 pop: 102372
},
{
 city: 'Erie',
 state: 'Pennsylvania',
 lat: 42.12922409,
 lon: -80.085059,
 pop: 101786
},

];

  var field = VectorField.read(windData, true);

var mapAnimator;
var legendSpeeds = [1, 3, 5, 10, 15, 30];

var MapMask = function(image, width, height) {
	this.image = image;
	this.width = width;
	this.height = height;
};

MapMask.prototype.endMove = function(animator) {
	this.move(animator);
}

MapMask.prototype.move = function(animator) {
	var s = this.image.style;
	s.width = ~~(animator.scale * this.width) + 'px';
	s.height = ~~(animator.scale * this.height) + 'px';
	s.left = animator.dx + 'px';
	s.top = animator.dy + 'px';
};

function isAnimating() {
	return document.getElementById('animating').checked;
}

function showCities() {
	document.getElementById('city-display').style.visibility =
	    document.getElementById('show-cities').checked ? 'visible' : 'hidden';
}

function doUnzoom() {
	mapAnimator.unzoom();
}

function format(x) {
	x = Math.round(x * 10) / 10;
	var a1 = ~~x;
	var a2 = (~~(x * 10)) % 10;
	return a1 + '.' + a2;	
}

function init() {
	loading = false;
	var timestamp = windData.timestamp || 'unknown on unknown';
	var parts = timestamp.split('on');
	var time = parts[0].trim();
	var day = parts[1].trim().replace(' 0', ' '); // correct "01" dates.;
	document.getElementById('update-time').innerHTML =
	   '<span id="day">' + day +'</span><br>' + time + ' EDT' +
		 '<br><span id="time-explanation">(time of forecast download)</span>';

	var avg = field.averageLength / 1.15; // knots --> miles per hour
  var max = field.maxLength / 1.15;
		document.getElementById('average-speed').innerHTML =
	    '<br>top speed: <b>' + format(max / 10) + ' m/s</b><br>' +
			'average: <b>' + format(avg / 10) + ' m/s</b>';

	var canvas = document.getElementById('display');
	var imageCanvas = document.getElementById('image-canvas');
	var mapProjection = new ScaledAlbers(1900., 0., canvas.height, -135.0, 30.0);

	var isMacFF = navigator.platform.indexOf('Mac') != -1 &&
	              navigator.userAgent.indexOf('Firefox') != -1;
	var isWinFF = navigator.platform.indexOf('Win') != -1 &&
	              navigator.userAgent.indexOf('Firefox') != -1;
	var isWinIE = navigator.platform.indexOf('Win') != -1 &&
	              navigator.userAgent.indexOf('MSIE') != -1;
	var numParticles = isMacFF || isWinIE ? 3500 : 5000; // slowwwww browsers
	var display = new MotionDisplay(canvas, imageCanvas, field,
	                                numParticles, mapProjection);

  // IE & FF Windows do weird stuff with very low alpha.
  if (isWinFF || isWinIE) {
		display.setAlpha(.05);
	}

  var navDiv = document.getElementById("city-display");
	var unzoom = document.getElementById('unzoom');
	mapAnimator = new Animator(navDiv, isAnimating, unzoom);
	mapAnimator.add(display);
	
	var mask = new MapMask(document.getElementById('mask'), 600, 700);
	mapAnimator.add(mask);
	
	var callout = document.getElementById('callout');
	var hovercard = new MotionDetails(navDiv, callout, field,
	                                  mapProjection, mapAnimator);

	var cityCanvas = document.getElementById('city-display');
	cityDisplay = new CityDisplay(cities, cityCanvas, mapProjection);
	mapAnimator.add(cityDisplay);
	cityDisplay.move();

  var legendAnimator = new Animator(null, isAnimating);

  // Scale for speed.
	// Numerator comes from size of map.
	// Denominator is knots vs. mph.
  var speedScaleFactor = 40 / 1.15;
	for (var i = 1; i <= legendSpeeds.length; i++) {
		var c = document.getElementById('legend' + i);
		var legendField = VectorField.constant(
		    legendSpeeds[i - 1] * speedScaleFactor, 0, 0, 0, c.width, c.height);
		var legend = new MotionDisplay(c, null, legendField, 30);
		// normalize so colors correspond to wind map's maximum length!
		legend.maxLength = field.maxLength * speedScaleFactor;
		legendAnimator.add(legend);
	}
	mapAnimator.start(40);
	legendAnimator.start(40);
}
