/* Расчетные функции. Сырые данные и диагностические признаки сохраняются отдельно от округленного результата. */
(function (root) {
  'use strict';
  const roundTo = (x, step) => Math.round((x + Number.EPSILON) / step) * step;
  const areaStep = x => x <= 10 ? .01 : x <= 20 ? .05 : x <= 100 ? .1 : x <= 200 ? .5 : 1;
  const area = (shape, a, b) => shape === 'round' ? Math.PI * a * a / 4 : a * b;
  function geometry(shape, sections) {
    if (!['round', 'flat'].includes(shape) || sections.length !== 3 || sections.some(s => !(s.a > 0) || (shape === 'flat' && !(s.b > 0)))) throw Error('Введите размеры в трёх сечениях.');
    const raw = sections.map(s => area(shape, s.a, s.b));
    const min = Math.min(...raw), mean = raw.reduce((a,b) => a+b, 0)/3;
    return {sections:raw, min, mean, f0:roundTo(min,areaStep(min)), fE:roundTo(mean,areaStep(mean))};
  }
  function regression(points) {
    const n=points.length;
    if (n<2) return null;
    const xbar=points.reduce((s,p)=>s+p.x,0)/n, ybar=points.reduce((s,p)=>s+p.y,0)/n;
    let xx=0,xy=0,yy=0;
    for (const p of points) {const dx=p.x-xbar,dy=p.y-ybar;xx+=dx*dx;xy+=dx*dy;yy+=dy*dy;}
    if (!(xx>0 && yy>0)) return null;
    const slope=xy/xx, intercept=ybar-slope*xbar, r2=xy*xy/(xx*yy);
    return {slope,intercept,r2,n};
  }
  function offsetYield(data, fit, offset=.002) {
    if (!fit || !(fit.slope>0)) return null;
    let previous=null;
    for (const p of data) {
      if (p.x < offset) continue;
      const difference=p.y-(fit.intercept+fit.slope*(p.x-offset));
      if (previous && previous.difference>0 && difference<=0) {
        const t=previous.difference/(previous.difference-difference);
        return {stress:previous.p.y+t*(p.y-previous.p.y),strain:previous.p.x+t*(p.x-previous.p.x)};
      }
      previous={p,difference};
    }
    return null;
  }
  // Автоматическое предложение оператору: спад усилия >= 0,5% до максимальной нагрузки.
  // Исключаем посадку образца и шум ниже 10% Pmax; подтверждение оператором обязательно.
  function physicalCandidate(data) {
    if (data.length<8) return null;
    const max=Math.max(...data.map(p=>p.force));
    const limit=data.findIndex(p=>p.force===max);
    for(let i=2;i<limit-2;i++) {
      if(data[i].force < .1*max || data[i].force < data[i-1].force || data[i].force < data[i-2].force) continue;
      let bottom=null;
      for(let j=i+1;j<Math.min(i+150,limit);j++) {
        if(!bottom || data[j].force<bottom.force) bottom=data[j];
        if(data[j].force>=data[i].force) break;
      }
      if(bottom && (data[i].force-bottom.force)/data[i].force>=.005 && bottom.x>data[i].x)
        return {upper:data[i].y,lower:bottom.y,upperStrain:data[i].x,lowerStrain:bottom.x};
    }
    return null;
  }
  function analyze(input) {
    const {shape,sections,l0,lk,neck,breakPosition,forceData,strainSource,gaugeBase,fitMin,fitMax,timeUnit}=input;
    const g=geometry(shape,sections);
    if(!(l0>0) || !Array.isArray(forceData) || forceData.length<10) throw Error('Нужны l₀ и минимум 10 точек диаграммы.');
    const warnings=[];
    const data=forceData.map((p,i)=>{
      if(!Number.isFinite(p.force) || p.force<0 || !Number.isFinite(p.deformation)) throw Error(`Некорректная точка ${i+1}.`);
      const x=strainSource==='percent' ? p.deformation/100 : strainSource==='ratio' ? p.deformation : p.deformation/gaugeBase;
      return {x,y:p.force/g.f0,force:p.force,time:p.time};
    });
    if(data.some((p,i)=>i && p.x<data[i-1].x-1e-9)) warnings.push('Деформация убывает: проверьте последовательность и прибор.');
    if(data.some((p,i)=>i && p.time!=null && data[i-1].time!=null && p.time<=data[i-1].time)) warnings.push('Время не возрастает строго; скорость не оценивается.');
    if(data.some(p=>!Number.isFinite(p.x))) throw Error('Некорректная деформация или база экстензометра.');
    const pmax=Math.max(...data.map(p=>p.force)), sigmaB=pmax/g.f0;
    let fit=null, e=null, offset=null;
    if(strainSource!=='crosshead' && Number.isFinite(fitMin) && Number.isFinite(fitMax) && fitMax>fitMin) {
      const subset=data.filter(p=>p.y>=fitMin && p.y<=fitMax && p.force>0);
      // Для E используется среднее значение площади; регрессия относится к тому же набору точек.
      const candidate=regression(subset.map(p=>({x:p.x,y:p.force/g.fE})));
      if(candidate && candidate.n>=50 && candidate.slope>0 && candidate.r2>=.9995) {
        e=candidate;
        const strengthFit=regression(subset.map(p=>({x:p.x,y:p.y})));
        offset=offsetYield(data,strengthFit);
        if(!offset) warnings.push('Пересечение линии смещения 0,2 % не найдено на записи.');
      } else warnings.push('E и σ₀.₂ не выданы: требуется ≥50 точек в выбранном упругом участке и r² ≥ 0,9995.');
    } else warnings.push('E и σ₀.₂ не определяются по ходу траверсы или без заданного диапазона регрессии.');
    let delta=null,psi=null;
    if(Number.isFinite(lk) && lk>=l0) {
      if(breakPosition==='near') { delta=100*(lk-l0)/l0; warnings.push('Разрыв в пределах 1/3 l₀ от метки: если δ ниже требования к продукции, примените перенос места разрыва по 7.6.3; прямой результат требует проверки.'); }
      else if(breakPosition==='ok') delta=100*(lk-l0)/l0;
      else warnings.push('Положение разрыва не подтверждено; δ не выдано.');
    } else warnings.push('Нет корректной конечной расчётной длины; δ не выдано.');
    if(neck && neck.a>0 && neck.b>0) {
      const fk=shape==='round'?Math.PI*neck.a*neck.b/4:neck.a*neck.b;
      if(fk<=g.f0) psi=100*(g.f0-fk)/g.f0;
      else warnings.push('Fк больше F₀; ψ не выдано.');
    } else warnings.push('Не введены два размера минимального сечения после разрыва; ψ не выдано.');
    const phys=physicalCandidate(data);
    if(phys) warnings.push('Физическая текучесть: найден кандидат по спаду усилия ≥0,5 %. Требуется проверка диаграммы оператором.');
    const hasTime=data.every(p=>Number.isFinite(p.time));
    let speed=null;
    if(hasTime && !warnings.some(w=>w.includes('Время не возрастает'))) {
      const subset=data.filter(p=>p.y>=fitMin&&p.y<=fitMax);
      const rate=regression(subset.map(p=>({x:p.time,y:p.x})));
      const load=regression(subset.map(p=>({x:p.time,y:p.y})));
      if(rate && load && subset.length>=2) speed={strainRate:rate.slope/(timeUnit==='min'?60:1),stressRate:load.slope/(timeUnit==='min'?60:1),points:subset.length};
    }
    if(!speed) warnings.push('Скорость не оценена: нужны время и точки в выбранном интервале.');
    return {g,data,pmax,sigmaB,e,offset,delta,psi,phys,speed,warnings};
  }
  const api={roundTo,areaStep,geometry,regression,offsetYield,physicalCandidate,analyze};
  root.GOST1497=api;
  if(typeof module!=='undefined') module.exports=api;
})(typeof window!=='undefined'?window:globalThis);
