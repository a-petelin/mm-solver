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
	function Solvers() {
	}	
	Solvers.EllipticSolver = function(mesh) {
	    
	}
	
	Solvers.gm = 5/3;
	Solvers.G = 6.67e-11;
	
	Solvers.SignalSpeed = function(xL, xR, i, Pst) {
		var S = {};
		//S["L"] = xL["v"+i] - Math.sqrt((Solvers.gm-1)*((xL.E+xL.p)/xL.ro-xL["v"+i]*xL["v"+i]/2));
		S["L"] = xL["v"+i] - xL["a"]*((Pst > xL["p"])?Math.sqrt(1+(1+Solvers.gm)/(2*Solvers.gm)*(Pst/xL["p"] - 1)):1);
		if(Math.abs(S["L"]) > 3e8) throw new Error("SL gone above light speed!");
		//S["R"] = xR["v"+i] + Math.sqrt((Solvers.gm-1)*((xR.E+xL.p)/xR.ro-xR["v"+i]*xR["v"+i]/2));
		S["R"] = xR["v"+i] + xR["a"]*((Pst > xR["p"])?Math.sqrt(1+(1+Solvers.gm)/(2*Solvers.gm)*(Pst/xR["p"] - 1)):1);
		if(Math.abs(S["R"]) > 3e8) throw new Error("SR gone above light speed!");
		return S;
	}

	Solvers.RiemanSolver = function(xL, xR, i, dimn) {
		var Pst = Math.max(0,(xL["p"]+xR["p"])/2-(xL["ro"]+xR["ro"])/2*(xL["a"]+xR["a"])/2*(xR["v"+i] - xL["v"+i])/2);
		var S = Solvers.SignalSpeed(xL, xR, i, Pst);	
		var SL = S.L;
		var SR = S.R;
		var Sst = (xR["p"] - xL["p"] + xL["ro"]*xL["v"+i]*(SL-xL["v"+i]) - xR["ro"]*xR["v"+i]*(SR-xR["v"+i]))/(xL["ro"]*(SL-xL["v"+i])-xR["ro"]*(SR-xR["v"+i]));
		var FLro = xL["ro"]*xL["v"+i];		
		var FRro = xR["ro"]*xR["v"+i];
		var ULro = xL["ro"];
		var URro = xR["ro"];
		var FLv = new Array(dimn);//xL["ro"]*xL["v"+i]*xL["v"+i] + xL["p"];
		var FRv = new Array(dimn);//xR["ro"]*xR["v"+i]*xR["v"+i] + xR["p"];
		var ULv = new Array(dimn);//xL["ro"]*xL["v"+i];
		var URv = new Array(dimn);//xR["ro"]*xR["v"+i];
		for (var k = 0 ; k < dimn; k++) {
			FLv[k] = xL["ro"]*xL["v"+i]*xL["v"+k];
			FRv[k] = xR["ro"]*xR["v"+i]*xR["v"+k];
			if(k == i) {
				FLv[k] += xL["p"];
				FRv[k] += xR["p"];
			}
			ULv[k] = xL["ro"]*xL["v"+k];
			URv[k] = xR["ro"]*xR["v"+k];
		}
		var FLE = xL["v"+i]*(xL["E"] + xL["p"]);
		var FRE = xR["v"+i]*(xR["E"] + xR["p"]);
		var ULE = xL["E"];
		var URE = xR["E"];
		var x = {};
		x.S = Math.max(Math.abs(SL), Math.abs(SR));
		if(SL > 0) {	
			x["Fro"] = FLro;
			for (var k = 0 ; k < dimn; k++) {
				x["Fv"+k] = FLv[k];
			}
			x["FE"] = FLE;
		} else if(SL <= 0 && Sst >= 0) {						
			var UstL = xL["ro"]*(SL - xL["v"+i])/(SL - Sst);
			x["Fro"] = FLro + SL*(UstL*1 - ULro);
			for (var k = 0 ; k < dimn; k++) {
				x["Fv"+k] = FLv[k] + SL*(UstL*(k==i?Sst:xL["v"+k]) - ULv[k]);
			}
			x["FE"] = FLE + SL*(UstL*(xL["E"]/xL["ro"]+(Sst-xL["v"+i])*(Sst+(xL["p"])/(xL["ro"]*(SL - xL["v"+i])))) - ULE);
		} else if(Sst <= 0 && SR >= 0) {						
			var UstR = xR["ro"]*(SR - xR["v"+i])/(SR - Sst);
			x["Fro"] = FRro + SR*(UstR*1 - URro);
			for (var k = 0 ; k < dimn; k++) {
				x["Fv"+k] = FRv[k] + SR*(UstR*(k==i?Sst:xR["v"+k]) - URv[k]);
			}
			x["FE"] = FRE + SR*(UstR*(xR["E"]/xR["ro"]+(Sst-xR["v"+i])*(Sst+(xR["p"])/(xR["ro"]*(SR - xR["v"+i])))) - URE);
		} else if(SR < 0) {
			x["Fro"] = FRro;
			for (var k = 0 ; k < dimn; k++) {
				x["Fv"+k] = FRv[k];
			}
			x["FE"] = FRE;
		}
		return x;
	}
	
	function Mesh(d, h, initf) {
		this.stepCounter = 0;
		this.isCriterialPoints = false;
		this.dims = d;
		this._h = h;
		this.blocks = [];		
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
		
		this._flowPrm = (function(){
			var initPrm = {ro: 0, v: 1, E: 0, p: 0, a: 0};
		    var flowParam = {};
			for(var prm in initPrm) {			
				if(initPrm[prm] == 1) {
					for(var i = 0; i < self.dims; i++){
						flowParam[prm+i] = "";
					}
				} else {
					flowParam[prm] = "";
				}
			}
		    return flowParam;
		})();
		
		this._data = {};
		
		this.dt = 1;
		
		if(initf !== undefined) {// STUB for test only
		    var c = [];
			for(var i=0;i<this.dims;i++) {
				c.push(0);
				this._data["_n-"+i] = [];
				this._data["_n"+i] = [];
			}
			for(var prm in this._movedParam) {
				this._data[prm] = [];
			}
			var dd = initf(c);
			for(var prm in this._calcPrm) {
				this._data[prm][0] = dd[prm];
			}	
			return;
		// Конструирование полной сетки
		    
		}
	}
	Mesh.prototype._moveNewToOld = function() {	
		var self = this;
		this.foreach(0,function(d,i){
			for(var prm in self._movedParam) {
				if((d[prm+"_new"]||[])[i] !== undefined) {
					d[prm][i] = d[prm+"_new"][i];
					delete d[prm+"_new"][i];
				}
			}
		});
	};
	Mesh.prototype.foreach = function(i,f, b, c) {	
		// STUB for test only
		f(this._data, 0);
		return;
		
	    b= b||this._f[i];
	    c = c||this._sizeData(this.dims);
	    for(var i = b, k = 0; i<this._sizeData(this.dims) && k < c; i = this._data["_i"+(n%this.dims)][i], k++){
		    f(this._data, i, this._data["_i"+(n%this.dims)][i]);
	    }
	    return i;
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
	function MMesh(d,k,h,curn, initf) {
		this.k = k;
		this.prevStep = STEP.FORWARD;
		this.i = 0;
		this.mm = [new Mesh(d,h, initf)];
		this.end = false;
		var self = this;
		this.kcr = 1/k/k;
		this.curn = curn;
		this.stepers = 
		{
				FORWARD: function (){
					var m = self._mesh();
					// 0. (new) E(new)?=>: (old)
					m._moveNewToOld();
					// а. (old) Построение решения гиперболических уравнений по вычисленным потокам (new)
					if(m.Smax && !self.i) {
					    m.dt = self.curn * m._h / m.Smax;
					}
					m.foreach(0, function(d,i){
					    var d_new = {};
					    var dfi = [];
					    for(var j = 0; j < m.dims; j++) {
    						var j_1 = d["_n-"+j][i];
    						var j1 = d["_n"+j][i];
    						j_1 = (j_1 === undefined?i:j_1);
    						j1 = (j1 === undefined?i:j1);
    						var r = (j_1==j?0:1)+(j1==j?0:1);
							if(r > 0) {
								dfi[j] = (((d.fi||[])[j1]||0)-((d.fi||[])[j_1]||0))/(r*m._h);
							} else {
								dfi[j] = 0;
							}
                        }
						for(var prm in m._calcPrm) {							
							d_new[prm] = 0;
							for(var j = 0; j < m.dims; j++) {
							    var f1,f_1;
							    f_1 = d["_F"+prm+"_"+j+"_1"]||[];
							    f1 = d["_F"+prm+"_"+j+"_2"]||[];
							    f_1 = f_1[i]||0;
							    f1 = f1[i]||0;
							    d_new[prm] += m.dt/m._h*(f_1 - f1);
							    switch(prm[0]) {
							    case "v":
							        d_new[prm] -= m.dt*d["ro"][i]*dfi[j];
							        break;
							    case "E":
							        d_new[prm] -= m.dt*d["ro"][i]*d["v"+j][i]*dfi[j];
							        break;
							    }
							}
						}
						if((d["ro"][i] + d_new["ro"]) <= 0) throw new Error("ro gone below zero!");	
				        for (var k = 0 ; k < m.dims; k++) {
					        if(Math.abs((d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][i] + d_new["ro"])) > 3e8) throw new Error("v"+k+" gone above light speed!");
							if(d["v"+k+"_new"] === undefined) d["v"+k+"_new"] = [];
					        d["v"+k+"_new"][i] = (d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][i] + d_new["ro"]);
				        }
				        if((d["E"][j] + d_new["E"]) <= 0) throw new Error("E gone below zero!");
						if(d["E_new"] === undefined) d["E_new"] = [];
				        d["E_new"][i] = d["E"][i] + d_new["E"];
						if(d["ro_new"] === undefined) d["ro_new"] = [];
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
		                if((Solvers.gm-1)*(d.E_new[i]-d.ro_new[i]*vv/2) < 0) throw new Error("p gone below zero!");
		                if(d["p_new"] === undefined) d["p_new"] = [];
				        d["p_new"][i] = (Solvers.gm-1)*(d.E_new[i]-d.ro_new[i]*vv/2);	
		                if(d["a_new"] === undefined) d["a_new"] = [];
				        d["a_new"][i] = Math.sqrt(Solvers.gm*d.p_new[i]/d.ro_new[i]);
					    
						d["_fcr"][i] = m._h/Math.sqrt(Math.PI*d.a_new[i]*d.a_new[i]/Solvers.G/d.ro_new[i]);
						
						var fn_cr = -Infinity;
						for(var prm in m._calcPrm) {
							prm += "_new";
							for(var j =0; j < m.dims; j++) {
								var j_1 = d["_n-"+j][i]
								var j1 = d["_n"+j][i];
								j_1 = (j_1 === undefined?i:j_1);
								j1 = (j1 === undefined?i:j1);
								var fnx_i = d[prm][j_1]+d[prm][j1];
								fn_cr = Math.max(fn_cr, (fnx_i==0)?Math.abs(d[prm][i]):Math.abs(2*d[prm][i]/(fnx_i)-1));
							}
						}
						d["_gcr"][i] = fn_cr;						
					});
					// в. Проставить сеточный флаг "есть точки по критериям"
					m.foreach(0, function(d,i){
						if(d["_fcr"][i] > 0.25 || d["_gcr"][i] > self.kcr) {
							var isMonade = (d["_is_monade"]||[])[i];
							if(!isMonade) {
								m.isCriterialPoints = true;
							}
							if(d["_is_move_down"] === undefined) d["_is_move_down"] = [];
							d["_is_move_down"][i] = true;
						}
					})
					// г. (new) Решение эллиптических уравнений (new)
					for(var i = 0; i < m.blocks.length; i++) {
					    Solvers.EllipticSolver(m.blocks[i]);
					}
					// д. (new) Рассчёт потоков для гиберболических уравнений (new)
					m.Smax = 0;
					for(var i = 0; i < m.dims; i++) {
					    m.foreach(i, (function(){var p_flow;
					            return function(d,j){
					                var j_1 = d["_n-"+i][j]
								    var j1 = d["_n"+i][j];
								    j_1 = (j_1 === undefined?j:j_1);
								    j1 = (j1 === undefined?j:j1);
								    var x={},x_1={},x1={};
								    for(var prm in m._flowPrm) {
								        x[prm] = d[prm+"_new"][j];
								        x_1[prm] = d[prm+"_new"][j_1];
								        x1[prm] = d[prm+"_new"][j1];
								    }
								    var flow = Solvers.RiemanSolver(x, x1, i, m.dims);
								    if(!p_flow) {
    							        p_flow = Solvers.RiemanSolver(x_1, x, i, m.dims);
    							    }
								    for(var k = 1; k < 3; k++) {
	    							    for(var prm in m._calcPrm) {
	    							        if(d["_F"+prm+"_"+i+"_"+k+"_new"] === undefined) d["_F"+prm+"_"+i+"_"+k+"_new"] = [];
	    							        d["_F"+prm+"_"+i+"_"+k+"_new"][j] = p_flow["F"+prm];
	    							    }
	    							    m.Smax = Math.max(m.Smax, p_flow.S);
	    							    p_flow = flow;
	    							}
					            };
					    })());
					}
					// е. Уменьшить счётчик шагов на 1
					m.stepCounter -= 1;
				},
				DOWN: function (){
					var m = self._mesh();
					// а. (old) Помещение точек согласно g,f-критериям в буфер, ели они есть
					var buff = [];
					m.foreach(0,function(d,i){
					    if(d["_is_move_down"][i]) {
					        var p = {};
							for(var prm in m._movedParam) {
							    if(d[prm][i] !== undefined) {
    								p[prm] = d[prm][i];
    							}
							}
							p.monade = i;
					        buff[i] = p;
					    }
					});
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
							p[prm] = 0;
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
	
	var mm = new MMesh(2, 2, 1e17, 0.495, function(c) {
	    var d = {};
	    var r2 = 0;
	    for(var i = 0; i < c.length; i++) {
	        d["v"+i] = 0;
	        r2 += c[i]*c[i];
	    }
	    d.ro = 1e-21 + (1e-23/(1+Math.pow(1e-14*Math.sqrt(r2),4)));
	    d.E = d.ro*10;
	    return d;
	});
	mm._mesh().stepCounter = 0;
	do {
		mm.doStep(mm.selectStep());
		mm.drawStep([{canvas: null, prm: "ro"}]);					
	} while(!mm.end);
	console.log(mm);
})();
