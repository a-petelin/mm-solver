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

(function(_export){
    var kk = 0.25;
    function Solvers() {
    }    
    Solvers.sineTransform = function(real) {
        // Generate complex vector with odd-symmetry
        var tmpImag = [0];
        var tmpReal = [0];
        var N = 2 * real.length + 2;
        for(var i = 0; i <= real.length; i++) {
            tmpReal[i+1] = real[i]||0;
            tmpReal[N-i-1] = -real[i]||0;
            tmpImag[i+1] = 0;
            tmpImag[N-i-1] = 0;
        }        
        // Perform FFT 
        transform(tmpReal, tmpImag);
        // Result of sine transform is imaginary part of the transformed vector
        for(var i = 0; i < real.length; i++) {
            real[i] = -Math.sqrt(1/N)*tmpImag[i+1];
        }
    } 
    
    Solvers.ExecutionError = function(message, step) {
        this.name = 'ExecutionError';
        this.message = message || 'Сообщение по умолчанию';
        this.step = step;
        this.stack = (new Error()).stack || "";
        this.stack = this.name +": " +this.message + "\n" + this.stack.split("\n").splice(2).join("\n");
    }
    Solvers.ExecutionError.prototype = Object.create(Error.prototype);
    Solvers.ExecutionError.prototype.constructor = Solvers.ExecutionError;
    
    Solvers.EllipticSolver = function(mesh, old) {
        var j = 0;
        var y = 0;
        var dat = new Array(mesh._N[0]);
        var W = [];
        var Wi = [];
        var Sk = [];
        var ki = [];
        for(var i = 0; i < mesh.dims; i++) {
            W[i] = Math.Complex(0,Math.PI/(2.0*(mesh._N[i]+1))).exp();
            Wi[i] = W[i];
            ki[i] = 1;
            Sk[i] = {"-1": Wi[i].re, "1": Wi[i].im};
        }
        mesh.foreach(0,function(d,i) {
            dat[j++] = i;
            if(j == mesh._N[0]) {
                var arrR = new Array(j), arrI = new Array(j);
                for (var k = 0; k < j; k++) { 
                    // Вместо самой функции u = 4*Pi*G*ro вычисляем невязку по приближённому значению fi включающему границы
                    arrR[k] = 4*Math.PI*Solvers.G*d.getDataValue("ro" + (!old?"_new":""),dat[k],0);
                    var LFi = -2*mesh.dims*d.getDataValue("_round_fi",dat[k],0);
                    var bFi = d.getDataValue("_border_fi",dat[k],[]);
                    var dFi = d.getAllNeib(dat[k]);
                    var test = 0;
                    for(var i = 0; i < 2*mesh.dims; i++) {
                        if(bFi[i] !== [][0]) {
                            LFi += bFi[i];
                        } else if(dFi[i] !== [][0]) {
                            LFi += d.getDataValue("_round_fi",dFi[i],0);
                        }
                    }
                    LFi /= mesh._h*mesh._h;
                    d.setDataValue("test",dat[k],arrR[k]-LFi);
                    arrR[k] -= LFi;
                }
                Solvers.sineTransform(arrR);
                for (var k = 0; k < j; k++) {
                    d.setDataValue("fst",dat[k],arrR[k]);
                }
                j = 0;
            }
        });
        for(var k = 1; k < mesh.dims ; k++) {
            j = 0;
            dat = new Array(mesh._N[k]);
            mesh.foreach(k, function(d, i) {
                dat[j++] = i;
                if(j == mesh._N[k]) {                        
                    var arrR = new Array(j), arrI = new Array(j);
                    for (var kk = 0; kk < j; kk++) {
                        arrR[kk] = d.getDataValue("fst",dat[kk]);
                    }
                    Solvers.sineTransform(arrR);
                    for (var kk = 0; kk < j; kk++) {
                        d.setDataValue("fst",dat[kk],arrR[kk]);
                    }
                    j = 0;
                }
            });
        }
        mesh.foreach(0, function(d, i, ni) {
            var denom = 0;
            for(var j = 0; j < mesh.dims; j++) {
                denom += -2*mesh.dims*(Sk[j][ki[j]]*Sk[j][ki[j]])/(mesh._h*mesh._h);
            }
            if(Math.abs(denom) > 1e-50) {
                d.setDataValue("fst",i,d.getDataValue("fst",i)/denom);
            } else {
                d.setDataValue("fst",i,0);
            }
            
            for(var j = 0; j < mesh.dims; j++) {
                if(d.getDataValue("_c"+j,ni) != d.getDataValue("_c"+j,i)) {
                    Wi[j] = Wi[j].mul(W[j]);
                    Sk[j] = {"-1": Wi[j].re, "1": Wi[j].im};
                    var pj = j - 1;
                    if(Wi[pj] != [][0]) {
                        Wi[pj] = Wi[pj].mul(W[pj]);
                        ki[pj] *= -1;
                        Sk[pj] = {"-1": Wi[pj].re, "1": Wi[pj].im};
                    }
                }
            }
        });
        for(var i = 0; i < mesh.dims ; i++) {
            j = 0;
            dat = new Array(mesh._N[i]);
            mesh.foreach(i, function(d, k) {
                dat[j++] = k;
                if(j == mesh._N[i]) {
                    var arrR = new Array(j), arrI = new Array(j);
                    for (var kk = 0; kk < j; kk++) {
                        arrR[kk] = d.getDataValue("fst",dat[kk]);
                    }
                    Solvers.sineTransform(arrR);
                    for (var kk = 0; kk < j; kk++) {
                        d.setDataValue("fst",dat[kk],arrR[kk]);
                    }
                    j = 0;
                }
            });
        }
        mesh.foreach(0, function(d, i) {
            // Т.к. вместо самой функции мы использовали невязку, то получили мы функцию ошибки и её надо прибавить к исходнику невязки
            var fi = d.getDataValue("_round_fi",i,0) + d.getDataValue("fst",i);
            d.setDataValue("fi_new",i, fi>0?0:fi);
        });
        mesh.foreach(0, function(d, i) {
            var func = 4*Math.PI*Solvers.G*d["ro" + (!old?"_new":"")][i];
            var LFi = -2*mesh.dims*d.getDataValue("fi_new",i,0);
            var bFi = d.getDataValue("_border_fi",i,[]);
            var dFi = d.getAllNeib(i);
            var test = 0;
            for(var j = 0; j < 2*mesh.dims; j++) {
                if(bFi[j] !== [][0]) {
                    LFi += bFi[j];
                } else if(dFi[j] !== [][0]) {
                    LFi += d.getDataValue("fi_new",dFi[j],0);
                }
            }
            LFi /= mesh._h*mesh._h;
            d.setDataValue("_func_fi",i,func);
            d.setDataValue("_L_fi",i,LFi);
            d.setDataValue("_r_fi",i,func - LFi);
        });
    }
    
    Solvers.gm = 5/3;
    Solvers.G = 6.67e-11;
    
    Solvers.SignalSpeed = function(xL, xR, i, Pst) {
        var S = {};
        //S["L"] = xL["v"+i] - Math.sqrt((Solvers.gm-1)*((xL.E+xL.p)/xL.ro-xL["v"+i]*xL["v"+i]/2));
        S["L"] = xL["v"+i] - xL["a"]*((Pst > xL["p"])?Math.sqrt(1+(1+Solvers.gm)/(2*Solvers.gm)*(Pst/xL["p"] - 1)):1);
        if(Math.abs(S["L"]) > 3e8) throw new Solvers.ExecutionError("SL gone above light speed!");
        //S["R"] = xR["v"+i] + Math.sqrt((Solvers.gm-1)*((xR.E+xL.p)/xR.ro-xR["v"+i]*xR["v"+i]/2));
        S["R"] = xR["v"+i] + xR["a"]*((Pst > xR["p"])?Math.sqrt(1+(1+Solvers.gm)/(2*Solvers.gm)*(Pst/xR["p"] - 1)):1);
        if(Math.abs(S["R"]) > 3e8) throw new Solvers.ExecutionError("SR gone above light speed!");
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
    
    function Mesh(d, h, n, kcr, initf) {
        this.stepCounter = 0;
        this.isCriterialPoints = false;
        this.dims = d;
        this._h = h;
        this.blocks = [];
        this.kcr = kcr;
        var self = this;

        this._movedParam = (function(){
            var initP = {v:2, ro:1, E:1, p:0, fi:0};
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
        this._critPrm = {ro:"", E:"", p:""};
        
        this._data = {
            setDataValue: function(prm, i, value) {
                if(this[prm] === undefined) this[prm] = [];
                this[prm][i] = value;
            },
            getDataValue: function(prm, i, undefValue) {
                var val = (this[prm]||[])[i];
                return val === [][0]?undefValue:val;
            },
            getNeibIndex: function(i, j, isReplaceI) {
                var j_1 = this["_n-"+j][i];
                var j1 = this["_n"+j][i];
                if(isReplaceI) {
                    j_1 = (j_1 === undefined?i:j_1);
                    j1 = (j1 === undefined?i:j1);
                }
                return [j_1,j1];
            },
            getAllNeib: function(i) {
                var nb = [];
                for(var j = 0; j < self.dims; j++) {
                    nb.push(this["_n-"+j][i]);
                    nb.push(this["_n"+j][i]);
                }
                return nb;
            }
        };
        for(var prm in this._movedParam) {
            this._data[prm] = [];
        }
        
        this._N = [];
        this._f = [];
        this.size = Math.pow(n, this.dims);
        
        if(initf !== undefined) {
            // Конструирование полной сетки
            var c = [], cc = [];
            for(var j = 0; j<this.dims; j++) {
                cc.push(0);
                this._N.push(n);
            }
            for(var i = 0; i < Math.pow(n, this.dims); i++) {
                for(var j = 0; j<this.dims; j++) {
                    c[j] = cc[j]*h-(n-1)*h/2;
                    this._data.setDataValue("_c"+j,i,c[j]);
                }
                this._data.setDataValue("fi",i,0);
                var dd = initf(c);
                for(var prm in this._calcPrm) {
                    this._data[prm][i] = dd[prm];
                }
                var d = 1, j = 0;
                do {
                    cc[j] += d;
                    d = 0;
                    if(cc[j] == n) {
                        cc[j] = 0;
                        d = 1;
                        j++;
                    }
                }while(d!=0);
            }
            this.composeBlocks();
            for(var i = 0; i < this.blocks.length; i++) {
                if(this.blocks[i]._data) {
                    Solvers.EllipticSolver(this.blocks[i], true);
                }
            }
        }
    }
    Mesh.prototype._moveNewToOld = function() {    
        var self = this;
        this.foreach(0,function(d,i){
            for(var prm in self._movedParam) {
                if((d[prm+"_new"]||[])[i] !== undefined) {
                    d.setDataValue(prm,i,d[prm+"_new"][i]);
                    delete d[prm+"_new"][i];
                }
            }
        });
        this.isCriterialPoints = false;
    };
    Mesh.prototype._dropNew = function() {
        var self = this;
        this.foreach(0,function(d,i){
            for(var prm in self._movedParam) {
                if((d[prm+"_new"]||[])[i] !== undefined) {
                    delete d[prm+"_new"][i];
                }
            }
        })
        this.isCriterialPoints = false;
    }
    Mesh.prototype.foreach = function(j,f) {
        for(var i = this._f[j]; i<this.size; i = this._data["_i"+(j%this.dims)][i]){
            f(this._data, i, this._data["_i"+(j%this.dims)][i]);
        }
    };
    Mesh.prototype.deletePoint = function(i) {        
        for(var p in this._data) {
            delete this._data[p][i];
        }
    }
    Mesh.prototype.addDataValue = function(prm, i, value) {
        this._data.setDataValue(prm,i,value);
    }
    Mesh.prototype.composeBlocks = function() {
        var jj = [], prev_i, cur_block_n;
        var delta = 1;
        var self = this;
        this.blocks = [];
        this._data[Object.keys(this._calcPrm)[0]].forEach(function(x,i){
            var ii = i, cj = [], ci = [];
            for(var j = 0; j < self.dims; j++) {
                ci[j] = ii % self._N[j];
                ii = (ii - ci[j]) / self._N[j];
            }
            for(var j = 1; j < self.dims; j++) {
                var arr = jj[j-1]||[];
                ci.push(ci.shift());
                var cj = 0;
                for(var k = self.dims-1; k > -1; k--) {
                    cj *= self._N[k];
                    cj += ci[k];
                }
                arr[cj] = i;
                jj[j-1] = arr;
            }
            if(prev_i === undefined) {
                self._f[0] = i;
                self.blocks.push(new MeshBlock(self));
                cur_block_n = 0;
                self._data.setDataValue("_block_b",i,true);
            } else {
                self._data.setDataValue("_i0",prev_i,i);
                var isNeib = (i - prev_i) == delta;
                for(var j = 1; j < self.dims; j++) {
                    isNeib = isNeib && (self._data["_c"+j][prev_i] == self._data["_c"+j][i]);
                }
                if(isNeib) {
                    self._data.setDataValue("_n0",prev_i,i);
                    self._data.setDataValue("_n-0",i,prev_i);
                } else {                            
                    self.blocks.push(new MeshBlock(self));
                    cur_block_n++;
                    self._data.setDataValue("_block_b",i,true);
                }
            }
            self.blocks[cur_block_n]._data.push(i);
            self.blocks[cur_block_n]._N[0] += 1;
            self.blocks[cur_block_n].size += 1;
            self._data.setDataValue("_block_n",i,cur_block_n);
            prev_i = i;
        });
        self._data.setDataValue("_i0",prev_i,self.size);
        for(var i = 1; i < self.dims; i++) {
            prev_i = [][0];
            delta *= self._N[i-1];
            jj[i-1].forEach(function(j,x){
                if(prev_i === undefined) {
                    self._f[i] = j;
                } else {
                    self._data.setDataValue("_i"+i,prev_i,j);
                    var isNeib = true;
                    for(var k = 0; k < self.dims; k++) {
                        if(k == i) {
                            isNeib = isNeib && (j - prev_i == delta);
                        } else {
                            isNeib = isNeib && (self._data["_c"+k][prev_i] == self._data["_c"+k][j]);
                        }
                    }
                    if(isNeib) {
                        self._data.setDataValue("_n"+i,prev_i,j);
                        self._data.setDataValue("_n-"+i,j,prev_i);
                        var block = self.blocks[self._data["_block_n"][j]];
                        if(self._data["_block_b"][j]){
                            var prev_block = self.blocks[self._data["_block_n"][prev_i]];
                            var isComp = true;
                            for(var k = 0; k < self.dims; k++) {
                                if(k != i) {
                                    isComp = isComp && block._N[k] == prev_block._N[k];
                                }
                            }
                            if(isComp) {
                                prev_block._N[i] += 1;
                                block = self._data["_block_n"][prev_i];
                                self.blocks[self._data["_block_n"][j]] = block;
                                self._data["_block_b"][j] = false;
                            }
                        }
                        if(!block._data) {
                            self._data["_block_n"][j] = block;
                            self.blocks[block]._data.push(j);
                            self.blocks[block].size += 1;
                        }
                    }
                }
                prev_i = j;
            });
            self._data.setDataValue("_i"+i,prev_i,self.size);
        }    
        for(var i = 0; i < self.blocks.length; i++) {
            if(self.blocks[i]._data) {   
                self.blocks[i]._data.sort(function(a, b) {
                    return a - b;
                });
                self._data["_block_b"][self.blocks[i]._data[0]] = false;
            }
        }
    }
    Mesh.prototype.setMoveDown = function(i) {
        var isMonade = (this._data["_is_monade"]||[])[i];
        if(!isMonade) {
            this.isCriterialPoints = true;
            this._data.setDataValue("_move_down",i,1);
        };
    }
    function MeshBlock(mesh) {
        this._m = mesh;
        this.dims = mesh.dims;
        this._N = ("0"+new Array(this.dims).join("1")).split("").map(function(x){return 1*x;});
        this._data = [];
        this.size = 0;
        this._h = mesh._h;
    }
    MeshBlock.prototype.foreach = function(i,f, b, c) {    
        var delta = 1;
        for(var j = 1; j <= i; j++) {
            delta *= this._N[j-1];
        }
        var next_j;
        for(var j = 0; j < this.size; j = next_j ) {
            next_j = (j + delta) % this.size + Math.floor((j + delta) / this.size);
            var next_i = this._data[next_j];
            if(j == this.size-1) {
                next_i = this._m.size;
                next_j = this.size;
            }
            f(this._m._data, this._data[j], next_i);
        }
    }
    MeshBlock.prototype.addDataValue = function(prm, i, value) {
        this._m._data.setDataValue(prm,i,value);
    }
    var STEP = {
                FORWARD: "FORWARD",
                DOWN: "DOWN",
                UP: "UP"
            };
    function MMesh(d,k,h,n, steps, curn, initf) {
        this.k = k;
        this.prevStep = STEP.FORWARD;
        this.i = 0;
        this.kcr = Math.pow(k, -2); 
        this.mm = [new Mesh(d,h, n, Math.pow(k, -2), initf)];
        this.mm[0].stepCounter = steps;
        this.mm[0].dt = 1;
        this.end = false;
        var self = this;
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
                            var r = (j_1==i?0:1)+(j1==i?0:1);
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
                                f1 =  d["_F"+prm+"_"+j+"_2"]||[];
                                f_1 = f_1[i]||0;
                                f1 =  f1[i]||0;
                                d_new[prm] += m.dt/m._h*(f_1 - f1);
                                switch(prm) {
                                case ("v"+j):
                                    d_new[prm] -= m.dt*d["ro"][i]*dfi[j];
                                    break;
                                case "E":
                                    d_new[prm] -= m.dt*d["ro"][i]*d["v"+j][i]*dfi[j];
                                    break;
                                }
                            }
                        }
                        if((d["ro"][i] + d_new["ro"]) <= 0) throw new Solvers.ExecutionError("ro gone below zero!", STEP.FORWARD);
                        for (var k = 0 ; k < m.dims; k++) {
                            if(Math.abs((d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][i] + d_new["ro"])) > 3e8) throw new Solvers.ExecutionError("v"+k+" gone above light speed!", STEP.FORWARD);
                            d.setDataValue("v"+k+"_new",i, (d["v"+k][i]*d["ro"][i] + d_new["v"+k])/(d["ro"][i] + d_new["ro"]));
                        }
                        if((d["E"][j] + d_new["E"]) <= 0) throw new Solvers.ExecutionError("E gone below zero!", STEP.FORWARD);
                        d.setDataValue("E_new",i,d["E"][i] + d_new["E"]);
                        d.setDataValue("ro_new",i,d["ro"][i] + d_new["ro"]);
                    });
                    // б. (new) Рассчёт g и f-критериев
                    m.foreach(0, function(d,i){
                        var vv = 0;
                        for(var j = 0; j < m.dims; j++) {
                            vv += d["v"+j+"_new"][i]*d["v"+j+"_new"][i];
                        }
                        if((Solvers.gm-1)*(d.E_new[i]-d.ro_new[i]*vv/2) < 0) throw new Solvers.ExecutionError("p gone below zero!", STEP.FORWARD);
                        d.setDataValue("p_new",i,(Solvers.gm-1)*(d.E_new[i]-d.ro_new[i]*vv/2));
                        d.setDataValue("a_new",i,Math.sqrt(Solvers.gm*d.p_new[i]/d.ro_new[i]));
                        
                        d.setDataValue("_fcr",i, m._h/Math.sqrt(Math.PI*d.a_new[i]*d.a_new[i]/Solvers.G/d.ro_new[i]));
                    });
                    m.foreach(0, function(d,i){
                        var fn_cr = -Infinity;
                        var L_cr = -Infinity;
                        for(var prm in m._critPrm) {
                            prm += "_new";
                            var _gr_eps = 1e-30;
                            for(var j =0; j < m.dims; j++) {
                                var j_1 = d["_n-"+j][i]
                                var j1 = d["_n"+j][i];
                                j_1 = (j_1 === undefined?i:j_1);
                                j1 = (j1 === undefined?i:j1);
                                var fnx_i = d[prm][j_1]+d[prm][j1];
                                fn_cr = Math.max(fn_cr, (Math.abs(fnx_i) < _gr_eps)?Math.abs(d[prm][i]):Math.abs(2*d[prm][i]/(fnx_i)-1));
                            }
                            
                            var L_fn = -2*m.dims*d[prm][i];
                            var g_fn_norm = 0;
                            var nbI = d.getAllNeib(i);
                            for(var j = 0; j < 2*m.dims; j++) {
                                L_fn += d.getDataValue(prm, nbI[j], d[prm][i]);
                                if(j % 2) {
                                    var fn1 = d.getDataValue(prm, nbI[j], d[prm][i]);
                                    var fn_1 = d.getDataValue(prm, nbI[j-1], d[prm][i]);
                                    var r = (nbI[j]==[][0]?0:1) + (nbI[j-1]==[][0]?0:1);
                                    g_fn_norm += Math.pow((fn1 - fn_1)/(r*m._h),2);
                                }
                            }            
                            L_fn = Math.abs(L_fn/m._h*m._h)/Math.pow(1+g_fn_norm, 3/2);
                            d.setDataValue("_L_"+prm,i,L_fn);
                        }
                        d.setDataValue("_gcr", i, fn_cr);
                    });
                    // в. Проставить сеточный флаг "есть точки по критериям"
                    //if(m.canDown) {
                        m.foreach(0, function(d,i){
                            var isBuff = (d["_buff_p"]||[])[i];
                            if(!isBuff && (d["_fcr"][i] > 0.125 || d["_gcr"][i]*0 > m.kcr)) {
                                m.setMoveDown(i);
                            }
                        });
                    //}
                    // г. (new) Решение эллиптических уравнений (new)
                    for(var i = 0; i < m.blocks.length; i++) {
                        if(m.blocks[i]._data) {
                            Solvers.EllipticSolver(m.blocks[i]);
                        }
                    }
                    // д. (new) Рассчёт потоков для гиберболических уравнений (new)
                    m.Smax = 0;
                    for(var i = 0; i < m.dims; i++) {
                        m.foreach(i, (function(){
                            var p_flow;
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
                                        d.setDataValue("_F"+prm+"_"+i+"_"+k+"_new",j,p_flow["F"+prm]);
                                    }
                                    m.Smax = Math.max(m.Smax, p_flow.S);
                                    p_flow = flow;
                                }
                                if(j1 == j) {
                                    p_flow = [][0];
                                }
                            };
                        })());
                    }
                    // е. Уменьшить счётчик шагов на 1
                    m.stepCounter -= 1;
                },
                DOWN: function (){
                    var m = self._mesh();
                    var buff = [];
                    buff.getPoint = function(m,d,i) {
                        var p = {};
                        for(var prm in m._movedParam) {
                            if(d[prm][i] !== undefined) {
                                p[prm] = d[prm][i];
                            }
                        }
                        for(var j=0;j<m.dims;j++) {
                            p["_c"+j] = d["_c"+j][i];
                        }
                        p._monade = i;
                        return p;
                    }
                    buff.addRealPoint = function(m,d,i) {
                        this[i] = this.getPoint(m,d,i);
                    }
                    buff.addBuffPoint = function(k,m,d,i) {
                        if(!d["_move_down"][i]) {
                            var p = this.getPoint(m,d,i);
                            p._buff_p = true;
                            this[i] = p;
                            d["_move_down"][i] = k;
                        }
                    }
                    // а. (old) Помещение точек согласно g,f-критериям в буфер, если они есть
                    m.foreach(0,function(d,i){
                        if(d["_move_down"][i]) {
                            buff.addRealPoint(m,d,i);
                            d.setDataValue("_is_monade",i,true);
                        } else if((d["_is_monade"]||[])[i]) {
                            d["_move_down"][i] = 1;
                        }
                    });
                    // б. (old) Помещение граничных точек в буфер (граничная точка сетки -- это такая точка не удолетворяющая g,f-критериям, но имеющая таковую в соседях по Муру на расстоянии 1)
                    for(var i=0;i<m.dims;i++) {
                        m.foreach(i,function(d,j){
                            if(d["_move_down"][j] == (i+1)) {
                                var j_1 = d["_n-"+i][j];
                                var j1 = d["_n"+i][j];
                                j_1 = (j_1 === undefined?j:j_1);
                                j1 = (j1 === undefined?j:j1);
                                buff.addBuffPoint(i+2,m,d,j_1);
                                buff.addBuffPoint(i+2,m,d,j1);
                                d["_move_down"][j] = i+2;
                            }
                        });
                    }
                    m.foreach(0,function(d,i){
                        d["_move_down"][i] = 0;
                    });
                    // в. Разбить каждую точку в буфере на куб из k^d точек.
                    buff.forEach(function(p,i, buff){
                        var cube = [];
                        var ii = [];
                        var jj = p._monade;
                        var kk = [];
                        var dd = [];
                        for(var l = 0; l < m.dims; l++) {
                            ii[l] = jj % m._N[l];
                            jj = (jj - ii[l]) / m._N[l];
                            ii[l] *= self.k;
                            kk.push(0);
                            dd.push((dd[l-1]*m._N[l-1]*self.k)||1);
                        }
                        for(var j = 0;j<Math.pow(self.k, m.dims);j++){
                            var p1 = {};
                            for(var prm in p) {
                                p1[prm] = p[prm];
                            }
                            p1.i = 0;
                            for(var k=0;k<m.dims;k++) {
                                p1.i += (ii[k]+kk[k])*dd[k];
                                p1["_c"+k] = (p["_c"+k] - 0.5*m._h) + (kk[k]+0.5)*m._h/self.k;
                            }
                            cube.push(p1);
                            var d = 1;
                            var e = 0;
                            do {
                                kk[e] += d;
                                d = 0;
                                if(kk[e] == self.k) {
                                    kk[e] = 0;
                                    d = 1;
                                    e++;
                                }
                            }while(d!=0);
                        }
                        buff[i] = cube;
                    });
                    // г. Увеличить номер текущей сетки на 1
                    var dt = m.dt / self.k;
                    self.i++;
                    m = self._mesh();
                    m.dt = dt;
                    // д. (old) Вставить точки из буфера в сетку
                    buff.forEach(function(cube,i, buff){
                        for(var j = 0; j < cube.length; j++) {
                            var p = cube[j];
                            var k = p.i;
                            delete p.i;
                            for(var prm in p) {
                                m._data.setDataValue(prm,k,p[prm]);
                            }
                            delete cube[j];
                        }
                        delete buff[i];
                    });
                    buff = {};
                    delete buff;
                    // е. (old) Разбить сетку на блоки
                    m.composeBlocks();
                    // ё. (old) Для всех точек каждого блока забрать по "монадам" решения эллиптических уравнений в качестве приближённых
                    var prev_m = self.mm[self.i-1];
                    m.foreach(0, function(d,i){
                        d.setDataValue("_round_fi", i, prev_m._data["fi"][d["_monade"][i]]);
                        for(var j = 0; j < m.dims; j++) {
                            var jj = [d["_n-"+j][i], d["_n"+j][i]];
                            for(var k = 0; k < 2; k++) {
                                if(jj[k] === undefined) {
                                    var ind = prev_m._data["_n"+(k > 0 ? "" : "-")+j][d["_monade"][i]];
                                    if(ind !== undefined) {
                                        var _border_fi = d.getDataValue("_border_fi",i,[]);
                                        _border_fi[j*2+k] = prev_m._data["fi"][ind];
                                        d.setDataValue("_border_fi",i,_border_fi);
                                    }
                                }
                            }
                        }
                    });
                    // *. 
                    for(var i = 0; i < m.blocks.length; i++) {
                        if(m.blocks[i]._data) {
                            Solvers.EllipticSolver(m.blocks[i],true);
                        }
                    }
                    // 0. (new) => (old)
                    prev_m._moveNewToOld();
                    // ж. Установить счётчик шагов на k
                    m.stepCounter = self.k;
                },
                UP:  function (){
                    var m = self._mesh();
                    // 0. (new) => (old)
                    m._moveNewToOld();
                    // а. (old) Поместить все неграничные точки сетки в буфер
                    var buff = [];
                    m.foreach(0,function(d,i){
                        if(!((d["_buff_p"]||[])[i])) {
                            var monade = d["_monade"][i];
                            var bu = buff[monade]||[];
                            var p = {};
                            for(var prm in m._movedParam) {
                                if(prm.substr(0,2) == "_F") {
                                    var dd = (prm.substr(-1)=="2"?"":"-")+prm.substr(-3,1);
                                    var ni = d["_n"+dd][i];
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
                    buff.forEach(function(cube,i, buff){
                        var p = {};
                        for(var prm in m._movedParam) {
                            var cnt = 0;
                            if(prm != "fi"){
                                p[prm] = 0;
                                for(var j = 0; j < cube.length; j++) {
                                    var val = cube[j][prm];
                                    if(val !== undefined) {
                                        p[prm] += val;
                                        cnt++;
                                    }
                                }
                                p[prm] /= cnt;
                            }
                        }
                        buff[i] = p;
                    });
                    // Уменьшить номер текущей сетки на 1.
                    //var S = m.Smax;
                    self.i--;
                    m = self._mesh();
                    //m.Smax = Math.max(m.Smax,S);
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
                                                if(prm.substr(0,2) == "_F" && prm.substr(-1) == (k+1) && prm.substr(-3,1) == j) {
                                                    var nd = prm.substr(-1) % 2 + 1;
                                                    var fprm = prm.substr(0,prm.length - 1)+nd;
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
        if(this.i == this.mm.length) {
            var m = this.mm[this.i-1];
            this.mm.push(new Mesh(m.dims,m._h/this.k));
            for(var i = 0; i < m.dims; i++) {
                this.mm[this.i]._N[i] = m._N[i]*this.k;
            }
            this.mm[this.i].size = m.size * Math.pow(this.k, m.dims);
            this.mm[this.i].kcr = m.kcr * this.k;
        }
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
            if(step == STEP.FORWARD && this.i == 0) {
                for(var i = 0; i < this.mm.length; i++) {
                    this.mm[i].canDown = true;
                }
            }
        }
    };
    
    _export.MMesh = MMesh;
    _export.ExecutionError = Solvers.ExecutionError;
})(this);
