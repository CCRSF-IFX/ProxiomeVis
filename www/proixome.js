(function() {
  function byIdOrSuffix(id) {
    return document.getElementById(id) || document.querySelector('[id$="' + id + '"]');
  }

  var rdsLoadElapsedTimer = null;
  var rdsLoadStartedAtMs = null;
  var rdsLoadElapsedPrefix = 'Elapsed';

  function formatRdsLoadElapsedTime(seconds) {
    seconds = Math.max(0, Math.floor(Number(seconds) || 0));
    var minutes = Math.floor(seconds / 60);
    var remainingSeconds = seconds % 60;
    var hours = Math.floor(minutes / 60);
    var remainingMinutes = minutes % 60;

    if (hours > 0) {
      return hours + 'h ' + String(remainingMinutes).padStart(2, '0') + 'm ' + String(remainingSeconds).padStart(2, '0') + 's';
    }
    if (minutes > 0) {
      return minutes + 'm ' + String(remainingSeconds).padStart(2, '0') + 's';
    }
    return seconds + 's';
  }

  function updateRdsLoadElapsedDisplay() {
    var elapsed = byIdOrSuffix('rds_load_elapsed');
    if (!elapsed || !rdsLoadStartedAtMs) {
      return;
    }
    var seconds = (Date.now() - rdsLoadStartedAtMs) / 1000;
    elapsed.textContent = rdsLoadElapsedPrefix + ': ' + formatRdsLoadElapsedTime(seconds);
  }

  function stopRdsLoadElapsedTimer() {
    if (rdsLoadElapsedTimer) {
      window.clearInterval(rdsLoadElapsedTimer);
    }
    rdsLoadElapsedTimer = null;
    rdsLoadStartedAtMs = null;
  }

  function startRdsLoadElapsedTimer(startedAtMs, prefix) {
    startedAtMs = Number(startedAtMs || 0);
    if (!isFinite(startedAtMs) || startedAtMs <= 0) {
      startedAtMs = rdsLoadStartedAtMs || Date.now();
    }

    rdsLoadStartedAtMs = startedAtMs;
    rdsLoadElapsedPrefix = prefix || 'Elapsed';
    if (rdsLoadElapsedTimer) {
      window.clearInterval(rdsLoadElapsedTimer);
    }
    updateRdsLoadElapsedDisplay();
    rdsLoadElapsedTimer = window.setInterval(updateRdsLoadElapsedDisplay, 1000);
  }

  function setRdsLoadState(message) {
    message = message || {};
    var state = message.state || 'idle';
    var loading = state === 'running' || message.disabled === true;
    var button = byIdOrSuffix('load_rds_path');
    var status = byIdOrSuffix('rds_load_status');
    var progress = byIdOrSuffix('rds_load_progress');
    var progressBar = byIdOrSuffix('rds_load_progress_bar');
    var elapsed = byIdOrSuffix('rds_load_elapsed');

    if (button) {
      button.disabled = loading;
      button.classList.toggle('disabled', loading);
      button.setAttribute('aria-disabled', loading ? 'true' : 'false');
      button.textContent = message.buttonLabel || (loading ? 'Loading...' : 'Load Data');
    }

    if (status) {
      status.className = 'rds-load-status ' + state;
      status.textContent = message.message || '';
    }

    if (progress) {
      progress.className = 'rds-load-progress ' + state;
      progress.setAttribute('aria-hidden', state === 'idle' ? 'true' : 'false');
    }

    if (progressBar) {
      var value = Number(message.progress || 0);
      if (!isFinite(value)) {
        value = 0;
      }
      value = Math.max(0, Math.min(100, value));
      progressBar.style.width = value + '%';
      progressBar.setAttribute('aria-valuenow', String(Math.round(value)));
    }

    if (loading) {
      startRdsLoadElapsedTimer(message.startedAtMs, message.elapsedPrefix || 'Elapsed');
    } else {
      stopRdsLoadElapsedTimer();
      if (elapsed) {
        elapsed.textContent = message.elapsedLabel || '';
      }
    }
  }

  function registerRdsLoadHandler() {
    if (!window.Shiny || !Shiny.addCustomMessageHandler || window.proxiomeRdsLoadHandlerRegistered) {
      return;
    }
    Shiny.addCustomMessageHandler('proxiome-rds-load-state', setRdsLoadState);
    window.proxiomeRdsLoadHandlerRegistered = true;
  }

  registerRdsLoadHandler();
  document.addEventListener('shiny:connected', registerRdsLoadHandler);
  document.addEventListener('click', function(event) {
    var target = event.target;
    if (!target || !target.closest) {
      return;
    }

    var button = target.closest('#load_rds_path, [id$="load_rds_path"]');
    if (!button || button.disabled) {
      return;
    }

    window.setTimeout(function() {
      var pathInput = byIdOrSuffix('rds_server_path');
      var label = 'RDS file';
      if (pathInput && pathInput.value) {
        var parts = pathInput.value.split(/[\\/]/);
        label = parts[parts.length - 1] || label;
      }
      setRdsLoadState({
        state: 'running',
        disabled: true,
        buttonLabel: 'Loading...',
        message: 'Starting RDS load for ' + label + '.',
        progress: 3,
        startedAtMs: Date.now(),
        elapsedPrefix: 'Elapsed',
        elapsedLabel: 'Elapsed: 0s'
      });
    }, 0);
  }, true);
})();

