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
	function Mesh(d) {
		this.stepCounter = 0;
		this.isCriterialPoints = false;
		this.dims = d;
		
		this._movedParam = (function(){
			var initP = {v:"", ro:"",E:"",_p:"",_fi:""};
			var movedParam = {};
			for(var prm in initP) {
				var p = {};
				p[prm] = "";
				switch(prm[0]) {
				case "v":
					p = {};
					for(var i = 0; i < this.dims; i++) {
						movedParam[prm + i] = ""; 
						p[prm + i] = ""; 
					}
					break;
				case "_":
					p = {};
				default:				
					movedParam[prm] = ""; 
					break;
				}
				for(var pp in p) {
					for(var i = 0; i < this.dims; i++) {
						movedParam["_F" + pp + "_" + i + "_1"] = "";
						movedParam["_F" + pp + "_" + i + "_2"] = ""; 
					}					
				}
			}
			return movedParam;
		})();
		
		this._calcPrm = [];
		
		this._data = {};
	}
	Mesh.prototype._moveNewToOld = function() {		
		m.foreach(0,function(d,i){
			for(var prm in this._movedParam) {
				if(d[prm+"_new"][i] !== undefined) {
					d[prm][i] = d[prm+"_new"][i];
					delete d[prm+"_new"][i];
				}
			}
		});
	};
	Mesh.prototype.foreach = function(i,f){
		
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
	function MMesh(d,k) {
		this.k = k;
		this.prevStep = STEP.FORWARD;
		this.i = 0;
		this.mm = [new Mesh(d)];
		this.end = false;
		this.stepers = 
		{
				FORWARD: function (){
					var m = this._mesh();
					// 0. (new) E(new)?=>: (old)
					m._moveNewToOld();
					// а. (old) Построение решения гиперболических уравнений по вычисленным потокам (new)
					m.foreach(0, function(d,i){
						for(var prm in m._calcPrm) {
							
						}
					});
					// б. (new) Рассчёт g и f-критериев 
					// в. Проставить сеточный флаг "есть точки по критериям"
					// г. (new) Решение эллиптических уравнений (new)
					// д. (new) Рассчёт потоков для гиберболических уравнений (new)
					// е. Уменьшить счётчик шагов на 1
					m.stepCounter -= 1;
				},
				DOWN: function (){
					var m = this._mesh();
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
					var m = this._mesh();
					// 0. (new) => (old)
					m._moveNewToOld();
					// а. (old) Поместить все неграничные точки сетки в буфер
					var buff = [];
					m.foreach(0,function(d,i){
						if(!d["_buff_p"][i]) {
							var monade = d["_monade"][i];
							var bu = buff[]||[];
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
							for(int j = 0; j < x.length(); j++) {
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
					this.i--;
					m = this._mesh();
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
	
	var mm = new MMesh(2);
	mm._mesh().stepCounter = 0;
	do {
		mm.doStep(mm.selectStep());
		mm.drawStep([{canvas: null, prm: "ro"}]);					
	} while(!mm.end);
})();