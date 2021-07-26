/*
    The "main.js" contains core JS component supporting the PERF-Analysis module.

    This web application is developed with Semantic-ui frame work.
    The charts are build using Apex-Chart js charting library.

    All the data for the report is derived from analyse_data.js
    data = { info: {SeqInfo: {}, RepInfo: {} } }

    PlotData is a dictionary with key as the repeat class and value as a dictionary
    PlotData: { REPEAT_CLASS: { LENGTH: FREQUENCY } }

*/

function ReverseString(str) {   
    str = str.split(''); str.reverse();
    return str.join('');
} 

function reverseComplement(motif){
    let complement = '';
    let basePairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'};
    for (let i=0; i<motif.length; i++) {
        complement += basePairs[motif[i]];
    }
    return ReverseString(complement);
}

const getCycles = function(motif) {
    cycles = []
    for (let i=0; i<motif.length; i++) {
        const newMotif = motif.slice(i, motif.length) + motif.slice(0, i);
        if (cycles.indexOf(newMotif) === -1) { cycles.push(newMotif); }
    }
    cycles.sort();
    return cycles
}

const buildCyclicalVariations = function(motif) {
    const cycles = getCycles(motif);
    const revCycles = getCycles(reverseComplement(motif));
    for (let i=0; i<revCycles.length; i++) {
        const revCycle = revCycles[i];
        if (cycles.indexOf(revCycle) === -1) { cycles.push(revCycle); }
    }
    return cycles;
}

const rangeLengths = function(min, max) {
    return [...Array(max-min+1).keys()].map(i => i + min);
}


const openOptions = function() {
    const DOM = document.getElementsByClassName('chart-options')[0];
    if (DOM.classList.contains('active')) { DOM.classList.remove('active'); DOM.classList.add('inactive'); }
    else { DOM.classList.remove('inactive'); DOM.classList.add('active'); }
}

const titleDropDown = function(elemClass) {
    const DOM = document.getElementsByClassName(elemClass)[0];
    if (DOM.classList.contains('active')) { DOM.classList.remove('active'); DOM.classList.add('inactive'); }
    else { DOM.classList.remove('inactive'); DOM.classList.add('active'); }
}

const selectDropDownOpen = function(elemClass) {
    const DOM = document.getElementsByClassName(elemClass)[0];
    if (DOM.classList.contains('active')) { DOM.classList.remove('active'); DOM.classList.add('inactive'); }
    else { DOM.classList.remove('inactive'); DOM.classList.add('active'); }
}


// Updating report data
for (const key in data.info.SeqInfo){$(`.value.${key}`).html(data.info.SeqInfo[key])};
for (const key in data.info.RepInfo){if (key != 'PlotData') { $(`.value.${key}`).html(data.info.RepInfo[key]); }};
// data options
const parameters = {
    repeatSelectionType: 'Order',
    orderedRepeatOrder: 'Top',
    orderedRepeatsOrderDataType: 'Frequency',
    orderedRepeatsNum: 10,
    selectedOrderedRepeats: [],
    selectedSelectionRepeats: [],
    repeatAbundanceDataType: 'Frequency',
    kmerGrouping: false,
    repeatKmerPercentage: false,
    showCyclicalVariations: false,
}

// repeatClassAbundanceData with Repeat class as the key and [bases, frequency] as the value.
const kmerNames = ["Monomer","Dimer","Trimer","Tetramer","Pentamer","Hexamer","Heptamer","Octamer","Nonamer","Decamer","Undecamer","Dodecamer","Tridecamer","Tetradecamer","Pentadecamer","Hexadecamer","Heptadecamer","Octadecamer","Nonadecamer","Icosamer","Uncosamer","Docosamer","Tricosamer","Tetracosamer","Pentacosamer","Hexacosamer","Heptacosamer","Octacosamer","Nonacosamer","Triacontamer","Untriacontamer","Dotriacontamer","Tritriacontamer","Tetratriacontamer","Pentatriacontamer","Hexatriacontamer","Heptatriacontamer","Octatriacontamer","Nonatriacontamer","Tetracontamer","Untetracontamer","Dotetracontamer","Tritetracontamer","Tetratetracontamer","Pentatetracontamer","Hexatetracontamer","Heptatetracontamer","Octatetracontamer","Nonatetracontamer","Pentacontamer"]
const plotData = data.info.RepInfo.PlotData;
const allRepClasses = data.info.RepInfo.AllRepClasses;
const repeatClassAbundanceData = {};
const cyclicalAbundanceData = {};
let selectedRepeats = [];
const chartTitles = ['REPEAT ABUNDANCE', 'K-MER DISTRIBUTION', 'LENGTH vs FREQUENCY'];
let current_chart = 0;

let globalMinLength = data.info.RepInfo.MinLength;
let globalRepeatMaxLength = 0;
allRepClasses.forEach( e => {
    const lengths = Object.keys(plotData[e]).map( d => parseInt(d));
    const maxLength = Math.max(...lengths);
    if (maxLength > globalRepeatMaxLength) { globalRepeatMaxLength = maxLength; }
    let frequency = 0; let bases = 0;
    lengths.forEach(l => { const f = plotData[e][l].reduce((a, b) => parseInt(a)+parseInt(b)); frequency += f; bases += f * l; })
    repeatClassAbundanceData[e] = [bases, frequency];
})
const sortedRepeatClassesByFrequency = data.info.RepInfo.AllRepClasses.map(x => x);
sortedRepeatClassesByFrequency.sort((a,b) => { return repeatClassAbundanceData[a][1] - repeatClassAbundanceData[b][1] });
const sortedRepeatClassesByBases = data.info.RepInfo.AllRepClasses.map(x => x);
sortedRepeatClassesByBases.sort((a,b) => { return repeatClassAbundanceData[a][0] - repeatClassAbundanceData[b][0] });


let allMotifLengths = allRepClasses.map(x => x.length); allMotifLengths.sort();
allMotifLengths = allMotifLengths.filter( (value, index, self) => { return self.indexOf(value) === index; } )
const allKmers = allMotifLengths.map(m => kmerNames[m-1]);
const allKmersAbundanceData = allKmers.map(x => [0, 0]);
allRepClasses.forEach(a => { 
    const kIndex = allKmers.indexOf(kmerNames[a.length-1]);
    allKmersAbundanceData[kIndex][0] += repeatClassAbundanceData[a][0]; 
    allKmersAbundanceData[kIndex][1] += repeatClassAbundanceData[a][1]; 
});

parameters.plotMinLength = globalMinLength;
parameters.plotMaxLength = 100;
let plotLengths = rangeLengths(parameters.plotMinLength, parameters.plotMaxLength);
document.getElementById('min-length').value = globalMinLength;
document.getElementById('min-length').max = globalRepeatMaxLength;
document.getElementById('min-length').min = globalMinLength;
document.getElementById('max-length').value = parameters.plotMaxLength;
document.getElementById('max-length').max = globalRepeatMaxLength;
document.getElementById('max-length').min = globalMinLength;

$("#repeat-select").multiSelect({
    selectableOptgroup: true,
    afterSelect: function(d){
        d.forEach(function(e){ if (parameters.selectedSelectionRepeats.indexOf(e) == -1) { parameters.selectedSelectionRepeats.push(e) } });
        parameters.repeatSelectionType='Select'
    },
    afterDeselect: function(d){
        d.forEach(element => { parameters.selectedSelectionRepeats.splice(parameters.selectedSelectionRepeats.indexOf(element), 1); }); 
        parameters.repeatSelectionType='Select'
    } 
});

const updateSelectedOrderRepeats = function() {
    if (parameters.orderedRepeatsOrderDataType === 'Frequency') { parameters.selectedOrderedRepeats = sortedRepeatClassesByFrequency.map(x => x); }
    else { parameters.selectedOrderedRepeats = sortedRepeatClassesByBases.map(x => x); }
    if ( parameters.orderedRepeatOrder == 'Top') { parameters.selectedOrderedRepeats.reverse(); }
    parameters.selectedOrderedRepeats = parameters.selectedOrderedRepeats.slice(0, parameters.orderedRepeatsNum);
}

const changeOrderRepeatOrder = function(order) {
    const DOM = document.getElementsByClassName('select-dropdown content repeats-order')[0];
    DOM.classList.remove('active'); DOM.classList.add('inactive');
    document.getElementsByClassName('select-dropdown title repeats-order')[0].innerHTML = order;
    parameters.orderedRepeatOrder = order;
    updateSelectedOrderRepeats(); parameters.repeatSelectionType='Order';
}

const orderedRepeatsOrderDataTypeChange = function(dataType) {
    const DOM = document.getElementsByClassName('select-dropdown content repeats-order-datatype')[0];
    DOM.classList.remove('active'); DOM.classList.add('inactive');
    document.getElementsByClassName('select-dropdown title repeats-order-datatype')[0].innerHTML = dataType;
    parameters.orderedRepeatsOrderDataType =  dataType;
    updateSelectedOrderRepeats(); parameters.repeatSelectionType = 'Order';
}

const orderedRepeatsNumChange = function() {
    parameters.orderedRepeatsNum = parseInt(document.getElementById('num-order-repeats').value);
    updateSelectedOrderRepeats(); parameters.repeatSelectionType='Order';
}

updateSelectedOrderRepeats();
let mainChart = echarts.init(document.getElementById('plot-area'));

const repeatAbundanceTitleSubText = function() {
    if (parameters.repeatSelectionType == 'Order') {
        return `${parameters.orderedRepeatOrder} ${parameters.orderedRepeatsNum} repeats by ${parameters.orderedRepeatsOrderDataType}\nData type: ${parameters.repeatAbundanceDataType} Length Range: ${parameters.plotMinLength}-${parameters.plotMaxLength}bp`
    } else {
        return `${parameters.selectedSelectionRepeats.length} selected repeats\nData type: ${parameters.repeatAbundanceDataType} Length Range: ${parameters.plotMinLength}-${parameters.plotMaxLength}bp`
    }
}

const repeatAbundanceBasicOptions  = {
    grid: { top: 80 },
    textStyle: { fontFamily: 'Ubuntu Mono'},
    title: { text: 'Repeat Abundance', left: 'center', textStyle: { fontSize: 20 }, subtext: repeatAbundanceTitleSubText(), subtextStyle: { fontSize: 14, color: 'black', lineHeight: 18 }},
    tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        backgroundColor: 'white', textStyle: { color: 'black', fontWeight: 'bold' },
    },
    xAxis: { data: [], axisTick: { show: true, interval: 0 }, axisLabel: { show: true, interval: 0, rotate: 45, textStyle: { fontSize: 16 } } },
    yAxis: { 
        nameLocation: 'middle',
        nameTextStyle: { fontWeight: 'bold', fontSize: 16 }, nameGap: 80, name: parameters.repeatAbundanceDataType,
        axisLabel: { show: true, textStyle: { fontSize: 16 }, }
    },
    series: [],
    toolbox: {
        right: 20,
        feature: {
          saveAsImage: {  title: 'Save as PNG', show: true,  type: 'png',  name: `${data.info.SeqInfo.FileName}_RepeatAbundance` },
          dataView: {  title: 'View data',  readOnly: true,  lang: ['Data View', 'Go  back'] },
        }
    },
}

const plotRepeatAbundance = function() {
    const repeatAbundanceOptions = {...repeatAbundanceBasicOptions};
    const dataIndex = (parameters.repeatAbundanceDataType === 'Frequency') ? 1 : 0;

    mainChart.clear();
    repeatAbundanceOptions.xAxis.data = selectedRepeats;
    repeatAbundanceOptions.yAxis.name = parameters.repeatAbundanceDataType;
    repeatAbundanceOptions.title.subtext = repeatAbundanceTitleSubText();
    
    const repClassLengthRangeAbundance = selectedRepeats.map( a =>  {
        let total = 0;
        const allLengths = Object.keys(plotData[a]).map(l => parseInt(l));
        allLengths.forEach( l => { 
            if ( l <= parameters.plotMaxLength && l >= parameters.plotMinLength) {
                if (dataIndex === 1) { total += plotData[a][l].reduce((a, b) => a+b)  }
                else { total += ((plotData[a][l].reduce((a, b) => a+b)) * l); }
            }
        })
        return total; 
    })

    repeatAbundanceOptions.series = [{ 
        type: 'bar',
        stack: 'Repeat Abundance',
        data: repClassLengthRangeAbundance, animation: false, 
        label: { 
                show: true,
                formatter: function(d) {
                    const percent = ((d.data/repeatClassAbundanceData[selectedRepeats[d.dataIndex]][dataIndex])*100).toFixed(2);
                    if (percent != 100) { return `${percent} %` }
                    else { return '' }
                }
            }
        },{
        type: 'bar',
        stack: 'Repeat Abundance',
        data: selectedRepeats.map( (a,i) => { return repeatClassAbundanceData[a][dataIndex] - repClassLengthRangeAbundance[i]; }),
        tooltip: { show: false },
        itemStyle: { opacity: 0.3 }, animation: false,
        label: { 
                show: true,
                formatter: function(d) {
                    if (d.data !== 0) {
                        const percent = ((d.data/repeatClassAbundanceData[selectedRepeats[d.dataIndex]][dataIndex])*100).toFixed(2);
                        return `${percent} %`;
                    } else { return ''; }
                }, color: 'black' 
            }
    }];
    mainChart.setOption(repeatAbundanceOptions);
}
selectedRepeats = parameters.selectedOrderedRepeats;
plotRepeatAbundance();

const lengthAbundanceTitleSubText = function() {
    if (parameters.repeatSelectionType == 'Order') {
        return `${parameters.orderedRepeatOrder} ${parameters.orderedRepeatsNum} repeats by ${parameters.orderedRepeatsOrderDataType}\nLength Range: ${parameters.plotMinLength}-${parameters.plotMaxLength}bp`
    } else {
        return `${parameters.selectedSelectionRepeats.length} selected repeats\nLength Range: ${parameters.plotMinLength}-${parameters.plotMaxLength}bp`
    }
}

const lengthAbundanceSeries = function() {
    return selectedRepeats.map( a => {
        return {
            type: 'line', name: a,
            data: plotLengths.map(l => {
                if (l in plotData[a]) { return plotData[a][l].reduce((a,b) => { return a+b; }) }
                else { return 0; }
            })
        }
    })
}

const lengthAbundanceBasicOptions  = {
    grid: { top: 120 },
    textStyle: { fontFamily: 'Ubuntu Mono'},
    title: { text: 'Length vs Abundance', left: 'center', textStyle: { fontSize: 20 }, subtext: lengthAbundanceTitleSubText(), subtextStyle: { fontSize: 14, color: 'black', lineHeight: 18 }},
    tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'line' },
        backgroundColor: 'white', textStyle: { color: 'black', fontWeight: 'bold' },
    },
    xAxis: { data: plotLengths, axisTick: { show: true, interval: 0 }, axisLabel: { show: true, rotate: 45, textStyle: { fontSize: 16 }, showMaxLabel: true } },
    yAxis: { 
        nameLocation: 'middle',
        nameTextStyle: { fontWeight: 'bold', fontSize: 16 }, nameGap: 80, name: parameters.repeatAbundanceDataType,
        axisLabel: { show: true, textStyle: { fontSize: 16 }, }
    },
    series: lengthAbundanceSeries(),
    legend: { top: 80, type: 'scroll', icon: 'circle', textStyle: { fontWeight: 'bold', fontSize: 14 } },
    toolbox: {
        right: 20,
        feature: {
          saveAsImage: {  title: 'Save as PNG', show: true,  type: 'png',  name: `${data.info.SeqInfo.FileName}_lengthAbundance` 
          },
          dataView: {  title: 'View data',  readOnly: true,  lang: ['Data View', 'Go  back'] 
          },
        }
    },
}

const plotLengthAbundance = function() {
    plotLengths = rangeLengths(parameters.plotMinLength, parameters.plotMaxLength);
    mainChart.clear();
    const lengthAbundanceOptions = {...lengthAbundanceBasicOptions};
    lengthAbundanceOptions.xAxis.data = plotLengths;
    lengthAbundanceOptions.yAxis.name = parameters.repeatAbundanceDataType;
    lengthAbundanceOptions.title.subtext = lengthAbundanceTitleSubText();
    lengthAbundanceOptions.series = lengthAbundanceSeries();
    mainChart.setOption(lengthAbundanceOptions);
}

const repeatDistributionTitleSubText = function() {
    if (parameters.repeatSelectionType == 'Order') {
        return `${parameters.orderedRepeatOrder} ${parameters.orderedRepeatsNum} repeats by ${parameters.orderedRepeatsOrderDataType}`
    } else {
        return `${parameters.selectedSelectionRepeats.length} selected repeats`
    }
}

const repeatDistributionSeries = function(kmers, kmersAbundance) {
    const series = [];
    const dataIndex = (parameters.repeatAbundanceDataType === 'Frequency') ? 1 : 0;
    selectedRepeats.forEach( a => {
        const values = new Array(kmers.length).fill(0);
        const kIndex = kmers.indexOf(kmerNames[a.length - 1]);
        const value = repeatClassAbundanceData[a][dataIndex];
        const percent = ((value/allKmersAbundanceData[kIndex][dataIndex])*100).toFixed(2)
        values[kIndex] = (parameters.repeatKmerPercentage) ? percent : value ;
        series.push({
            name: a, data: values, type: 'bar', stack: 'Kmer', 
            xAxisIndex: 0, yAxisIndex: 0, label: {
                show: true, position: 'inside',
                formatter: function(d) {
                    const di = (parameters.repeatAbundanceDataType === 'Frequency') ? 1 : 0;
                    if (d.data !== 0) { return d.seriesName; } else { return ''; }
                }
            }
        });
    })
    const total = (dataIndex === 0) ? parseInt(data.info.RepInfo.TotalRepBases) : parseInt(data.info.RepInfo.TotalRepFreq);
    kmers.forEach((k,i) => { 
        const value = (parameters.repeatKmerPercentage) ? ((kmersAbundance[i]/total)*100).toFixed(2) : kmersAbundance[i];
        series.push({
            type: 'bar', name: k, stack: 'Distribution', data:[value],
            label: {
                show: true, position: 'inside',
                formatter: function(d) { 
                    if (parameters.repeatKmerPercentage) {
                        return d.seriesName;
                    } else {
                        return `${d.seriesName}\n${((d.data/kmersAbundance.reduce((a, b) => { return a+b; }))*100).toFixed(2)} %`; 
                    }
                }
            },
            xAxisIndex: 1, yAxisIndex: 1
        }) 
    })
    
    if (parameters.repeatKmerPercentage) {
        kmers.forEach( (kmer, i) => {
            const remValues = kmers.map(k => 0);
            remValues[i] = (((allKmersAbundanceData[i][dataIndex]-kmersAbundance[i])/allKmersAbundanceData[i][dataIndex])*100).toFixed(2)
            series.push({
                name: kmer, data: remValues, type: 'bar', stack: 'Kmer', 
                xAxisIndex: 0, yAxisIndex: 0, label: {
                    show: true, position: 'inside',
                    formatter: function(d) {
                        if (d.data !== 0) { return `${d.data} %`; }
                        else { return ''; }
                    }, textStyle: { color: 'black' }
                }, tooltip: {show: false}, itemStyle: { color: '#303030', opacity: 0.3 }
            });
        })
        series.push({
            type: 'bar', name: 'rem', stack: 'Distribution', data:[(((total-(kmersAbundance.reduce((a,b) => { return a+b; })))/total)*100).toFixed(2)],
            label: {
                show: true, position: 'inside',
                formatter: function(d) { return `${d.data} %` }
            },
            xAxisIndex: 1, yAxisIndex: 1, tooltip: { show: false }, itemStyle: { color: '#303030', opacity: 0.3 }
        })
    }
    return series;
}

const repeatDistributionBasicOptions  = {
    textStyle: { fontFamily: 'Ubuntu Mono'},
    grid: [{ top: 80, width: '60%' }, { top: 80, width: '10%', right: 100 }],
    title: { text: 'Repeat Distribution', left: 'center', textStyle: { fontSize: 20 }, subtext: repeatDistributionTitleSubText(), subtextStyle: { fontSize: 14, color: 'black', lineHeight: 18 }},
    tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        backgroundColor: 'white', textStyle: { color: 'black', fontWeight: 'bold' },
        formatter: function(d) {
            if (d[0].axisIndex === 0) {
                const kmer = d[0].axisValue;
                let content = '';
                d.forEach(val => { 
                    if (kmerNames[val.seriesName.length-1] === kmer) { 
                        content += `
                        <div style="display: flex; flex-direction: row; flex-wrap: nowrap; justify-content: flex-start; margin: 3px 0px;">
                            <div style="background-color: ${val.color}; width: 16px; height: 16px; border-radius: 8px; margin-right: 3px; display: inline-block"></div>
                            <div style="display: inline-block">${val.seriesName}: ${val.data}</div>
                        </div>`;
                    }
                })
                return content;
            }
            else { 
                let content = '';
                d.forEach(val => { 
                    content += `
                    <div style="display: flex; flex-direction: row; flex-wrap: nowrap; justify-content: flex-start; margin: 3px 0px;">
                        <div style="background-color: ${val.color}; width: 16px; height: 16px; border-radius: 8px; margin-right: 3px; display: inline-block"></div>
                        <div style="display: inline-block">${val.seriesName}: ${val.data}</div>
                    </div>`;
                })
                return content;
            }
        }
    },
    xAxis: [
        { gridIndex: 0, data: [], axisTick: { show: true, interval: 0 }, axisLabel: { show: true, rotate: 45, textStyle: { fontSize: 16 }, showMaxLabel: true } },
        { gridIndex: 1, data: ['Total'], axisLabel: { show: true, textStyle: { fontSize: 16 } } },
    ],
    yAxis: [{ 
        gridIndex: 0, nameLocation: 'middle',
        nameTextStyle: { fontWeight: 'bold', fontSize: 16 }, nameGap: 80, name: parameters.repeatAbundanceDataType,
        axisLabel: { show: true, textStyle: { fontSize: 16 }, }
    },{ 
        gridIndex: 1, nameLocation: 'middle', 
        nameTextStyle: { fontWeight: 'bold', fontSize: 16 }, nameGap: 80, name: parameters.repeatAbundanceDataType,
        axisLabel: { show: true, textStyle: { fontSize: 16 }, }, position: 'right',
    }],
    legend: { show: false, type: 'scroll', icon: 'circle', top: 60, textStyle: { fontSize: 14, fontWeight: 'bold' } },
    toolbox: {
        right: 20,
        feature: {
          saveAsImage: {  title: 'Save as PNG', show: true,  type: 'png',  name: `${data.info.SeqInfo.FileName}_lengthAbundance` 
          },
          dataView: {  title: 'View data',  readOnly: true,  lang: ['Data View', 'Go  back'] 
          },
        }
    },
}

const plotRepeatDistribution = function() {
    mainChart.clear();
    const dataIndex = (parameters.repeatAbundanceDataType === 'Frequency') ? 1 : 0;
    let motifLengths = selectedRepeats.map(x => x.length); motifLengths.sort();
    motifLengths = motifLengths.filter( (value, index, self) => { return self.indexOf(value) === index; } )
    const kmers = motifLengths.map(m => kmerNames[m-1]);
    
    let kmersAbundance = new Array(kmers.length).fill(0);
    selectedRepeats.forEach(a => { kmersAbundance[kmers.indexOf(kmerNames[a.length-1])] += repeatClassAbundanceData[a][dataIndex]; })
    
    const repeatDistributionOptions = {...repeatDistributionBasicOptions};
    repeatDistributionOptions.xAxis[0].data = kmers;
    if (parameters.repeatKmerPercentage) { repeatDistributionOptions.yAxis[0].max = 100; repeatDistributionOptions.yAxis[1].max = 100 }
    else { repeatDistributionOptions.yAxis[0].max = Math.max(...kmersAbundance); repeatDistributionOptions.yAxis[1].max = kmersAbundance.reduce((a, b) => { return a+b; }); }
    repeatDistributionOptions.yAxis[0].name = (parameters.repeatKmerPercentage) ? 'Percentage' : parameters.repeatAbundanceDataType;
    repeatDistributionOptions.yAxis[1].name = (parameters.repeatKmerPercentage) ? 'Percentage' : parameters.repeatAbundanceDataType;
    repeatDistributionOptions.title.subtext = repeatDistributionTitleSubText();
    repeatDistributionOptions.series = repeatDistributionSeries(kmers, kmersAbundance);
    mainChart.setOption(repeatDistributionOptions);
}


const updateChart = function() {
    console.log('Plotting Chart!');
    if (parameters.repeatSelectionType === 'Order') { selectedRepeats = parameters.selectedOrderedRepeats; }
    else { selectedRepeats = parameters.selectedSelectionRepeats }
    if (current_chart === 0) {
        plotRepeatAbundance();
    } else if (current_chart === 1) {
        plotRepeatDistribution();        
    } else {
        plotLengthAbundance();
    }
    const DOM = document.getElementsByClassName('chart-options')[0];
    DOM.classList.remove('active'); DOM.classList.add('inactive');
}

const repeatAbundanceDataTypeChange = function(dataType) { parameters.repeatAbundanceDataType = dataType; updateChart(); }

const kmerGroupingChange = function() {  
    parameters.kmerGrouping = document.getElementById('kmer-group-checkbox').checked; 
    updateChart(); 
}

const showActualRepeatsChange = function() { 
    showActualRepeats = document.getElementById('show-actual-repeats').checked;
}

const updateLengthChange = function(type) {
    if (type === 'min') { parameters.plotMinLength = parseInt(document.getElementById('min-length').value); }
    else { parameters.plotMaxLength = parseInt(document.getElementById('max-length').value); }
    updateChart();
}

const showPercentagesChange = function() {
    parameters.repeatKmerPercentage = document.getElementById('percentage-checkbox').checked; 
    updateChart();
}

const chartNav = function(change, current) {
    if (current === 3) {
        current_chart += change;
        current_chart = (current_chart < 0) ? 3+current_chart : current_chart;
        current_chart = (current_chart > 2) ? 0 : current_chart;
    } else if (change === 3) {
        current_chart = current;
        const DOM = document.getElementsByClassName('title-dropdown chart content')[0];
        DOM.classList.remove('active'); DOM.classList.add('inactive');
    }
    document.getElementById('charts-title').innerHTML = chartTitles[current_chart];
    for (let a of document.getElementsByClassName('parameters')) { a.style.display = 'initial'; }
    if (current_chart === 0) {
        for (let a of document.getElementsByClassName('parameters')) { if (!(a.classList.contains('repeat-abundance'))) { a.style.display = 'none'; } }
    } else if (current_chart === 1) {
        for (let a of document.getElementsByClassName('parameters')) { if (!(a.classList.contains('repeat-distribution'))) { a.style.display = 'none'; } }
    } else {
        for (let a of document.getElementsByClassName('parameters')) { if (!(a.classList.contains('length-abundance'))) { a.style.display = 'none'; } }
    }
    updateChart();
}


// Static Chart of Kmer Distribution
const kmerChart = echarts.init(document.getElementById('kmer-chart'));
let allkmerSeries = [];
allKmers.forEach((k,i) => { allkmerSeries.push({
    type: 'bar', name: k, stack: 'Distribution', data:[allKmersAbundanceData[i][1]],
    label: {show: true, position: 'inside', formatter: function(d) { return `${((d.data/data.info.RepInfo.TotalRepFreq)*100).toFixed(2)} %`; }}
}) })
const kmerChartOptions = {
    textStyle: { fontFamily: 'Ubuntu Mono',},
    grid: { containLabel: true, bottom: 20, top: 40, left: 60, right: 60 },
    title: { text: 'K-mer Abundance Distribution', left: 'center' },
    tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
    },
    yAxis: { data: ['Distribution'], axisTick: { show: true, interval: 0 }, axisLabel: { show: true, interval: 0, rotate: 0 } },
    xAxis: { max: data.info.RepInfo.TotalRepFreq, min: 0, axisTick: { show: true }, axisLabel: { show: true } },
    series: allkmerSeries,
    toolbox: {
        right: 20,
        feature: {
          saveAsImage: {  title: 'Save as PNG', show: true,  type: 'png',  name: `${data.info.SeqInfo.FileName}_KmerDistribution` 
          },
          dataView: {  title: 'View data',  readOnly: true,  lang: ['Data View', 'Go  back'] 
          },
        }
    },
}
kmerChart.setOption(kmerChartOptions);
window.onresize = function() { kmerChart.resize(); mainChart.resize(); }