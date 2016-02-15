/*
1. Прямо
	0. (new) E(new)?=>: (old)
	а. (old) Построение решения гиперболических уравнений по вычисленным потокам (new)
	б. (new) Рассчёт g и f-критериев 
	в. Проставить сеточный флаг "есть точки по критериям"
	г. (new) Решение эллиптических уравнений (new)
	д. (new) Рассчёт потоков для гиберболических уравнений (new)
	е. Уменьшить счётчик шагов на 1
	ё. (new) Отрисовать сетку
2. Вниз
	а. (old) Помещение точек согласно g,f-критериям в буфер, ели они есть
	б. (old) Помещение граничных точек в буфер (граничная точка сетки -- это такая точка не удолетворяющая g,f-критериям, но имеющая таковую в соседях по Муру на расстоянии 1)
	*. При помещении точек вместо решений эллиптических уравнений забирается усреднение этого решения по граням.
	в. Разбить каждую точку в буфере на куб из k^d точек.
	0. (new) => (old)
	г. Увеличить номер текущей сетки на 1
	д. (old) Вставить точки из буфера в сетку
	е. (old) Разбить сетку на блоки
	ё. (old) Для точек каждого блока остаить только те значения эллиптических уравнений, которые лежать на граничных для блока гранях.
	ж. Установить счётчик шагов на k
	з. (old) Отрисовать сетку
3. Вверх
	0. (new) => (old)
	а. (old) Поместить все неграничные точки сетки в буфер
	б. Каждый куб из k^d точек буфера собрать в одну точку
	в. Удалить точки согласно g,f-критериям 	
	г. Уменьшить номер текущей сетки на 1.
	д. (old) Вставить точки из буфера
	е. (old) В точках граничных с вставленными, потоки в эти точки установить согласно значениям на соответствующих гранях вставленных точек
	ё. (old) Отрисовать сетку
	
Изначально считается что был переход прямо
Если предыдущий переход был прямо, то 
-- Если "есть точки по критериям", то делаем переход вниз
-- Иначе если существует сетка уровнем ниже, то -- переход вниз
-- Иначе если существует сетка уровнем выше и счётчик шагов равен 0, то -- переход вверх
-- Иначе если счетчик шагов равен 0, то завершение
-- Иначе -- переход прямо
Если предыдущий переход был вниз, то делаем переход прямо
Если предыдущий переход был вверх, то
-- Если существует сетка уровнем выше и счётчик шагов равен 0, то делаем переход вверх
-- Иначе если счетчик шагов равен 0, то завершение
-- Иначе -- переход прямо
*/

(function(){
	var kk = 0.25;
	function Mesh(d, h, initf) {
		this.stepCounter = 0;
		this.isCriterialPoints = false;
		this.dims = d;
		this._h = h;		
		var self = this;

		this._movedParam = (function(){
			var initP = {v:2, ro:1, E:1, p:0, fi:0, _fcr:0, _gcr:0};
			var movedParam = {};
			for(var prm in initP) {
				var p = {};
				p[prm] = "";
				var val = initP[prm];
				switch(val) {
				case 2:
					p = {};
					for(var i = 0; i < self.dims; i++) {
						movedParam[prm + i] = ""; 
						p[prm + i] = ""; 
					}
					break;
				case 0:
					p = {};
				default:				
					movedParam[prm] = ""; 
					break;
				}
				for(var pp in p) {
					for(var i = 0; i < self.dims; i++) {
						movedParam["_F" + pp + "_" + i + "_1"] = "";
						movedParam["_F" + pp + "_" + i + "_2"] = ""; 
					}					
				}
			}
			return movedParam;
		})();
		
		this._calcPrm = (function(){
			var initPrm = {ro: 0, v: 1, E: 0};
			var calcParam = {};
			for(var prm in initPrm) {			
				if(initPrm[prm] == 1) {
					for(var i = 0; i < self.dims; i++){
						calcParam[prm+i] = "";
					}
				} else {
					calcParam[prm] = "";
				}
			}
			return calcParam;
		})();
		
		this._data = {};
		
		if(initf !== undefined) {
		// Конструирование полной сетки
		    
		}
	}
	Mesh.prototype._moveNewToOld = function() {		
		this.foreach(0,function(d,i){
			for(var prm in this._movedParam) {
				if(d[prm+"_new"][i] !== undefined) {
					d[prm][i] = d[prm+"_new"][i];
					delete d[prm+"_new"][i];
				}
			}
		});
	};
	Mesh.prototype.foreach = function(i,f, b, c) {
	    /*b= b||this._f[i];
	    c = c||this._sizeData(this.dims);
	    for(var i = b, k = 0; i<this._sizeData(this.dims) && k < c; i = this._data["_i"+(n%this.dims)][i], k++){
		    f(this._data, i, this._data["_i"+(n%this.dims)][i]);
	    }
	    return i;*/
	};
	Mesh.prototype.deletePiont = function(i) {		
		for(var p in this._data) {
			delete this._data[p];
		}
	}
	var STEP = {
				FORWARD: "FORWARD",
				DOWN: "DOWN",
				UP: "UP"
			};
	function MMesh(d,k,h, initf) {
		this.k = k;
		this.prevStep = STEP.FORWARD;
		this.i = 0;
		this.mm = [new Mesh(d,h, initf)];
		this.end = false;
		var self = this;
		this.stepers = 
		{
				FORWARD: function (){
					var m = self._mesh();
					// 0. (new) E(new)?=>: (old)
					m._moveNewToOld();
					// а. (old) Построение решения гиперболических уравнений по вычисленным потокам (new)
					m.foreach(0, function(d,i){
					    var d_new = {};
					    var dfi = [];
					    for(var j = 0; j < m.dims; j++) {
    						var j_1 = d["_n-"+j][i]
    						var j1 = d["_n"+j][i];
    						j_1 = (j_1 === undefined?j:j_1);
    						j1 = (j1 === undefined?j:j1);
    						r = (j_1==j?0:1)+(j1==j?0:1);
    						dfi[i] = (d._fi[j1]-d._fi[j_1])/(r*mesh._h);
                        }
						for(var prm in m._calcPrm) {							
							d_new[prm] = 0;
							for(var j = 0; j < m.dims; j++) {
							    d_new[prm] += m.dt/m._h*(d["_F"+prm+"_"+j+"_1"][i] - d["_F"+prm+"_"+j+"_1"][i]);
							    switch(prm[0]) {
							    case "v":
							        d_new[prm] -= dt*d["ro"][i]*dfi[i];
							        break;
							    case "E":
							        d_new[prm] -= dt*d["ro"][i]*d["v"+j][i]*dfi[i];
							        break;
							    }
							}
						}
						if((d["ro"][i] + d_new["ro"]) <= 0) throw new Error("ro gone below zero!");	
				        for (var k = 0 ; k < m.dims; k++) {
					        if(Math.abs((d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][j] + d_new["ro"])) > 3e8) throw new Error("v"+i+" gone above light speed!");
					        d["v"+k+"_new"][i] = (d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][j] + d_new["ro"]);
				        }
				        if((d["E"][j] + d_new["E"]) <= 0) throw new Error("E gone below zero!");
				        d["E_new"][i] = d["E"][i] + d_new["E"];
				        d["ro_new"][i] = d["ro"][i] + d_new["ro"];
					});
					// б. (new) Рассчёт g и f-критериев
					m.foreach(0, function(d,i){
		                var vv = 0;
		                for(var j = 0; j < m.dims; j++) {
			                vv += d["v"+j+"_new"][i]*d["v"+j+"_new"][i];
		                }
		                /*
		                var cv = d.v0[j]/Math.sqrt(vv);
		                var sv = d.v1[j]/Math.sqrt(vv);
		                d.vv[j] = (vv==0)?(-2*Math.PI):(Math.atan(sv/cv)+((cv<0)?(Math.sign(sv)*Math.PI):0));
		                */
		                var gm = 5/3;
		                if((gm-1)*(d.E[i]-d.ro[i]*vv/2) < 0) throw new Error("p gone below zero!");
		                d["p_new"][i] = (gm-1)*(d.E_new[i]-d.ro_new[i]*vv/2);	
		                d["a_new"][i] = Math.sqrt(gm*d.p_new[i]/d.ro_new[i]);
					    
					});
					// в. Проставить сеточный флаг "есть точки по критериям"
					// г. (new) Решение эллиптических уравнений (new)
					// д. (new) Рассчёт потоков для гиберболических уравнений (new)
					// е. Уменьшить счётчик шагов на 1
					m.stepCounter -= 1;
				},
				DOWN: function (){
					var m = self._mesh();
					// а. (old) Помещение точек согласно g,f-критериям в буфер, ели они есть
					// б. (old) Помещение граничных точек в буфер (граничная точка сетки -- это такая точка не удолетворяющая g,f-критериям, но имеющая таковую в соседях по Муру на расстоянии 1)
					// *. При помещении точек вместо решений эллиптических уравнений забирается усреднение этого решения по граням.
					// в. Разбить каждую точку в буфере на куб из k^d точек.
					// 0. (new) => (old)
					// г. Увеличить номер текущей сетки на 1
					// д. (old) Вставить точки из буфера в сетку
					// е. (old) Разбить сетку на блоки
					// ё. (old) Для точек каждого блока остаить только те значения эллиптических уравнений, которые лежать на граничных для блока гранях.
					// ж. Установить счётчик шагов на k
				},
				UP:  function (){
					var m = self._mesh();
					// 0. (new) => (old)
					m._moveNewToOld();
					// а. (old) Поместить все неграничные точки сетки в буфер
					var buff = [];
					m.foreach(0,function(d,i){
						if(!d["_buff_p"][i]) {
							var monade = d["_monade"][i];
							var bu = buff[monade]||[];
							var p = {};
							for(var prm in m._movedParam) {
								if(prm.substr(2) == "_F") {
									var d = (prm.substr(-1)=="2"?"":"-")+prm.substr(-3,1);
									var ni = d["_n"+d][i];
									if(ni === undefined || d["_monade"][ni] != monade) {
										p[prm] = d[prm][i];
									}
								} else {
									p[prm] = d[prm][i];
								}
							}
							bu.push(p);
							buff[monade] = bu;
						}
					});
					// Каждый куб из k^d точек буфера собрать в одну точку
					buff.forEach(function(x,i, buff){
						var p = {};
						for(var prm in m._movedParam) {
							var cnt = 0;
							for(var j = 0; j < x.length(); j++) {
								var val = x[j][prm];
								if(val !== undefined) {
									p[prm] += val;
									cnt++;
								}
							}
							p[prm] /= cnt;
						}
						buff[i] = p;
					});
					// Уменьшить номер текущей сетки на 1.
					self.i--;
					m = self._mesh();
					// д. (old) Вставить точки из буфера
					m.foreach(0, function(d,i) {
						if(buff[i]) {
							for(var prm in buff[i]) {
								d[prm][i] = buff[i][prm];
							}
							delete buff[i];
						}
					});
					delete buff;
					// е. (old) В точках граничных с вставленными, потоки в эти точки установить согласно значениям на соответствующих гранях вставленных точек
					m.foreach(0, function(d,i) {
						if(!d["_is_monade"][i]) {
							for(var j = 0; j < m.dims; j++) {
								var nb = [d["_n-"+j][i], d["_n"+j][i]];
								for(var k = 0; k < 2; k++) {
									if(nb[k] !==undefined) {
										if(d["_is_monade"][nb[k]]) {
											for(var prm in d) {
												if(prm.substr(2) == "_F" && prm.substr(-1) == k+1) {
													var nd = prm.substr(-1) % 2 + 1;
													var fprm = prm.substr(0,prm.length() - 1)+nd;
													d[prm][i] = d[fprm][nb[k]];
												}
											}
										}
									}
								}
							}
						}
					});
					// в. Удалить точки согласно g,f-критериям
					m.foreach(0,function(d,i){
						// ?????
					});
				}
		};
	}
	MMesh.prototype._mesh = function () {
		return this.mm[this.i];
	};
	MMesh.prototype._iUpMesh = function() {
		return this.mm[this.i-1] !== undefined;
	};
	MMesh.prototype._iDownMesh = function() {
		return this.mm[this.i+1] !== undefined;
	};
	MMesh.prototype.selectStep = function(){
		var step;
		switch(this.prevStep) {
			case STEP.FORWARD:
				if(this._mesh().isCriterialPoints) {
					step = STEP.DOWN;
				} else if(this._iDownMesh()) {
					step = STEP.DOWN;
				} else if(this._iUpMesh() && this._mesh().stepCounter == 0) {
					step = STEP.UP;
				} else if(this._mesh().stepCounter == 0) {
					this.end = true;
				} else {
					step = STEP.FORWARD;
				}
				break;
			case STEP.DOWN:
				step = STEP.FORWARD;
				break;
			case STEP.UP:
				if(this._iUpMesh() && this._mesh().stepCounter == 0) {
					step = STEP.UP;
				} else if(this._mesh().stepCounter == 0) {
					this.end = true;
				} else {
					step = STEP.FORWARD;
				}
				break;
		}
		return step;
	};
	MMesh.prototype.doStep = function(step) {
		if(!this.end) {
			this.stepers[step]();
			this.prevStep = step;
		}
	};
	MMesh.prototype.drawStep = function(drawObjs) {
		for(var drObj in drawObjs) {
			switch(this.prevStep) {
				case STEP.FORWARD:
				// ё. (new) Отрисовать сетку
					break;
				case STEP.DOWN:
				case STEP.UP:
				// з. (old) Отрисовать сетку
				// ё. (old) Отрисовать сетку
					break;
			}
		}
	};
	
	var mm = new MMesh(2, 2, 1e17, function(c) {
	    var d = {};
	    var r2 = 0;
	    for(var i = 0; i < c.length; i++) {
	        d["v"+i] = 0;
	        r2 += c[i]*c[i];
	    }
	    d.ro = 1e-21 + (0.01/(1+1e-14*r2*r2));
	    d.E = d.ro*10;
	    return d;
	});
	mm._mesh().stepCounter = 1;
	do {
		mm.doStep(mm.selectStep());
		mm.drawStep([{canvas: null, prm: "ro"}]);					
	} while(!mm.end);
})();
