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


function shuffle_ids_and_cats(ids, cats) {
    ids_fin = [];
    cat_fin = [];
    idxs = Array(ids.length).fill().map((element, index) => index);
    idxs_shuffled = jsPsych.randomization.sampleWithoutReplacement(idxs, ids.length);
    for (var i = 0; i < idxs.length; i++) {
        ids_fin[i] = ids[idxs_shuffled[i]];
        cat_fin[i] = cats[idxs_shuffled[i]];
    }
    both = { ids: ids_fin, cat: cat_fin };
    return (both);
}


