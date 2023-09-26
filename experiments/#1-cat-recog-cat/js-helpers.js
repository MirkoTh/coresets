// get subject ID
function getQueryVariable(variable) {
    var query = window.location.search.substring(1);
    var vars = query.split("&");
    for (var i = 0; i < vars.length; i++) {
        var pair = vars[i].split("=");
        if (pair[0] == variable) { return pair[1]; }
    }
    return (false);
}


function saveData(filedata, task) {
    var filename = "./data/" + task + "-participant-" + participant_id + ".json";
    $.post("save_data.php", { postresult: filedata + "\n", postfile: filename })
}

async function saveSeveralData(filedata, task) {
    var filename = "./data/" + task + "-participant-" + participant_id + ".json";
    var n_data = filedata.length;
    for (var i = 0; i < n_data; i++) {
        $.post("save_data.php", { postresult: JSON.stringify(filedata[i]) + "\n", postfile: filename })
    }
}


function download(content, fileName, contentType) {
    var a = document.createElement("a");
    var file = new Blob([content], { type: contentType });
    a.href = URL.createObjectURL(file);
    a.download = fileName;
    a.click();
}

function prepare_recall(data_recall) {
    var trial_id_recall = data_recall.select("trial_id_recall");
    var set_size = data_recall.select("set_size");
    var stimuli = data_recall.select("stimuli");
    var responses = data_recall.select("recall");
    var n_correct = data_recall.select("accuracy");
    var rt = data_recall.select("rt");
    var data_recall_clean = {
        participant_id: participant_id,
        trial_id_recall: trial_id_recall,
        set_size: set_size,
        stimuli: stimuli,
        response: responses,
        n_correct: n_correct,
        rt: rt
    };
    return (data_recall_clean)
}

function prepare_processing(data_processing) {
    var trial_id_recall = data_processing.select("trial_id_recall");
    var trial_id_processing = data_processing.select("trial_id_processing");
    var set_size = data_processing.select("set_size");
    var accuracy = data_processing.select("accuracy");
    var rt = data_processing.select("rt");
    var data_clean = {
        trial_id_recall: trial_id_recall,
        trial_id_processing: trial_id_processing,
        set_size: set_size,
        accuracy: accuracy,
        rt: rt
    };
    return (data_clean)
}

function make_stimuli(n_trials_train, n_trials_recog, n_trials_generalize, which_ids_0, which_ids_1) {
    var ids_0_train = [];
    var ids_1_train = [];
    var ids_0_lures = [];
    var ids_1_lures = [];
    var ids_0_generalize = [];
    var ids_1_generalize = [];

    var cat_0_train = [];
    var cat_1_train = [];
    var cat_0_lures = [];
    var cat_1_lures = [];
    var cat_0_generalize = [];
    var cat_1_generalize = [];


    // training stimuli
    for (var i = 0; i < n_trials_train / 2; i++) {
        ids_0_train[i] = which_ids_0[i];
        ids_1_train[i] = which_ids_1[i];
        cat_0_train[i] = 0;
        cat_1_train[i] = 1;
    }

    // recognition lures
    var counter_recog = 0;
    for (var i = n_trials_train / 2; i < n_trials_recog / 2; i++) {
        ids_0_lures[counter_recog] = which_ids_0[i];
        ids_1_lures[counter_recog] = which_ids_1[i];
        cat_0_lures[counter_recog] = 0;
        cat_1_lures[counter_recog] = 1;
        counter_recog += 1;
    }

    // generalization stimuli
    var counter_generalize = 0;
    for (var i = n_trials_recog / 2; i < ((n_trials_recog + n_trials_generalize) / 2); i++) {
        ids_0_generalize[counter_generalize] = which_ids_0[i];
        ids_1_generalize[counter_generalize] = which_ids_1[i];
        cat_0_generalize[counter_generalize] = 0;
        cat_1_generalize[counter_generalize] = 1;
        counter_generalize += 1;
    }

    stimuli = {
        ids_train: ids_0_train.concat(ids_1_train),
        cat_train: cat_0_train.concat(cat_1_train),

        ids_lures: ids_0_lures.concat(ids_1_lures),
        cat_lures: cat_0_lures.concat(cat_1_lures),

        ids_generalize: ids_0_generalize.concat(ids_1_generalize),
        cat_generalize: cat_0_generalize.concat(cat_1_generalize),
    }


    return (stimuli);
}


function shuffle_items(ids, cats, is_old = null) {
    ids_fin = [];
    cat_fin = [];
    is_old_fin = [];
    idxs = Array(ids.length).fill().map((element, index) => index);
    idxs_shuffled = jsPsych.randomization.sampleWithoutReplacement(idxs, ids.length);
    for (var i = 0; i < idxs.length; i++) {
        ids_fin[i] = ids[idxs_shuffled[i]];
        cat_fin[i] = cats[idxs_shuffled[i]];
        if (is_old != null) {
            is_old_fin[i] = is_old[idxs_shuffled[i]];
        }
    }
    both = { ids: ids_fin, cat: cat_fin, is_old: is_old_fin };
    return (both);
}


function condition_and_ncategories() {

    document.getElementById("condition_id").innerHTML = condition_id
    document.getElementById("n_categories").innerHTML = n_categories
    clickStart('page0', 'page1')
}

function direct_to_study() {
    window.location.href = "study.html";
}

var flag = 0;
var instcounter = 0;
function instructioncheck(pg, pg_prev) {
    var ch1 = 0;
    var ch2 = 0;
    var ch3 = 0;
    var ch4 = 0;
    var ch5 = 0;
    //check if correct answers are provided
    if (document.getElementById('icheck1').checked) { var ch1 = 1; color('q1icheck1', 'green') }
    else { colorWrongAnswer("q1", 'red') }
    if (document.getElementById('icheck2').checked) { var ch2 = 1; color('q2icheck2', 'green') }
    else { colorWrongAnswer("q2", 'red') }
    if (document.getElementById('icheck3').checked) { var ch3 = 1; color('q3icheck3', 'green') }
    else { colorWrongAnswer("q3", 'red') }
    if (document.getElementById('icheck4').checked) { var ch4 = 1; color('q4icheck4', 'green') }
    else { colorWrongAnswer("q4", 'red') }
    if (document.getElementById('icheck5').checked) { var ch5 = 1; color('q5icheck5', 'green') }
    else { colorWrongAnswer("q5", 'red') }

    var checksum = ch1 + ch2 + ch3 + ch4 + ch5;
    var criterion = 5;

    // indicate correct answers
    ++flag;
    clickStart(pg, pg);
    change("check", "Continue")

    // page transition 
    if ((checksum === criterion) && (flag == 2)) {
        //if correct, continue 
        //begintrial();
        direct_to_study();
        // alert
        alert('Great, you have answered all of the questions correctly. The study will now start.');
    }
    else {
        if (flag == 2) {
            instcounter++;
            colorWrongAnswer("q1", '#333333')
            colorWrongAnswer("q2", '#333333')
            colorWrongAnswer("q3", '#333333')
            colorWrongAnswer("q4", '#333333')
            colorWrongAnswer("q5", '#333333')
            //if one or more answers are wrong, raise alert box
            alert('You have answered some of the questions wrong. Please try again.');
            // go back to instructions
            clickStart(pg, pg_prev);
            flag = 0;

        }
    }

}

//function to hide one html div and show another
function clickStart(hide, show) {
    document.getElementById(hide).style.display = 'none';
    document.getElementById(show).style.display = 'block';
    window.scrollTo(0, 0);
}


// color text
function color(id, col) {
    document.getElementById(id).style.color = col;
}

//changes inner HTML of div with ID=x to y
function change(x, y) {
    document.getElementById(x).innerHTML = y;
}
