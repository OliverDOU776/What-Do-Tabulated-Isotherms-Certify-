'use strict';

// No remote services or analytics are used. Videos replay recorded outputs.
for (const button of document.querySelectorAll('[data-copy]')) {
  button.addEventListener('click', async () => {
    const target = document.getElementById(button.dataset.copy);
    const status = document.getElementById('copy-status');
    try {
      await navigator.clipboard.writeText(target.textContent);
      button.textContent = 'Copied';
      status.textContent = 'Copied to clipboard.';
      setTimeout(() => { button.textContent = 'Copy'; }, 1800);
    } catch (_) {
      const selection = window.getSelection();
      const range = document.createRange();
      range.selectNodeContents(target);
      selection.removeAllRanges();
      selection.addRange(range);
      status.textContent = 'Text selected. Use your browser copy command.';
    }
  });
}

// The opening cartoon loops silently; visitors who prefer reduced motion get a paused video with controls.
const teaser = document.getElementById('teaser-video');
if (teaser && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
  teaser.removeAttribute('autoplay'); teaser.pause(); teaser.controls = true;
}
