/*
    The "tables.js" contains JS code which updates the data for
    tables present in the PERF-analysis module.

    This web application is developed with Semantic-ui frame work.
    The charts are build using Apex-Chart js charting library.

    All the data for the report is derived from analyse_data.js
    data = {info: {genomeInfo: {}, repInfo: {}, plotInfo: {}}}

*/


const table_types = ['REPEATS SUMMARY', 'LONGEST REPEATS', 'MOST UNITS'];
let current_table = 0;

const summaryTableData = [];

const sortedRepeats = data.info.RepInfo.AllRepClasses.map(x => x);
sortedRepeats.sort((a,b) => { return repeatClassAbundanceData[b][2] - repeatClassAbundanceData[a][2] });
console.log(sortedRepeats)
for (let rep in sortedRepeats) {
    rep = sortedRepeats[rep];
    summaryTableData.push({
        'repClass': rep,
        'freq': repeatClassAbundanceData[rep][1],
        'freqPercent': ((repeatClassAbundanceData[rep][1]/data.info.RepInfo.TotalRepeats)*100).toFixed(2),
        'bases': repeatClassAbundanceData[rep][0],
        'basesPercent': ((repeatClassAbundanceData[rep][0]/data.info.RepInfo.RepeatBases)*100).toFixed(2),
        'reads': repeatClassAbundanceData[rep][2],
        'readsPercent': ((repeatClassAbundanceData[rep][0]/data.info.RepInfo.RepeatReads)*100).toFixed(2)
    })
}

const summaryColumns = [ 'Repeat Class', 'Reads', '% Reads', 'Frequency', '% Frequency', 'Bases', '% Bases'];
const repeatLocColumns = ['Sequence id', 'Start', 'Stop', 'Repeat Class', 'Repeat length', 'Strand', 'Units', 'Actual Repeat'];
const summaryColKeys = ["repClass", "reads", "readsPercent", "freq", "freqPercent", "bases", "basesPercent"];
const repeatColKeys = ["seq", "start", "end", "repClass", "repLength", "repOri", "repUnit", "actualRep"];
const updateTableData = function(tableId, tableData, type) {
    let header = [];
    let colKeys = [];
    if (type === 'summary') { header = summaryColumns; colKeys = summaryColKeys; }
    else { header = repeatLocColumns; colKeys = repeatColKeys; }

    const tableDOM = document.getElementById(tableId);
    const table = document.createElement('table');
    // table.className = "ui sortable celled table";
    const tableHead = document.createElement('thead')
    const tableHeadRow = document.createElement('tr');
    header.forEach(function(e){ const headCell = document.createElement('th'); headCell.innerHTML = e; tableHeadRow.appendChild(headCell); })
    tableHead.appendChild(tableHeadRow);
    
    const tableBody = document.createElement('tbody');
    for (let d in tableData) {
        d = tableData[d];
        const row = document.createElement('tr');
        colKeys.forEach(function(e){ const cell = document.createElement('td'); cell.innerHTML = d[e]; row.appendChild(cell); })
        tableBody.appendChild(row);
    }
    table.appendChild(tableHead);
    table.appendChild(tableBody);
    tableDOM.appendChild(table);
}

const tableNav = function(change, current) {
    if (current === 3) {
        current_table += change;
        current_table = (current_table < 0) ? 3+current_table : current_table;
        current_table = (current_table > 2) ? 0 : current_table;
    } else if (change === 3) {
        current_table = current;
        const DOM = document.getElementsByClassName('title-dropdown table content')[0];
        DOM.classList.remove('active'); DOM.classList.add('inactive');
    }
    console.log(current_table)
    document.getElementById('data-table').innerHTML = '';
    document.getElementById('table-title').innerHTML = table_types[current_table];
    switch (current_table) {
        case 0: updateTableData('data-table', summaryTableData, 'summary');
        case 1: updateTableData('data-table', data.info.RepInfo.LongestRepeats, 'longest');
        case 2: updateTableData('data-table', data.info.RepInfo.MostRepeatUnits, 'units');
    }
}
updateTableData('data-table', summaryTableData, 'summary');
