/* Phase 13: Guided Tour Engine — TUTR-01..TUTR-05
   UI-SPEC §State Machine, §Six-Step Tour Config, §Copywriting Contract
   State: localStorage molaop_tour_seen / molaop_tour_step
   Cross-page: hidden #enrichmentForm tour field → server tour_active →
               <meta name="tour-active"> → resume at step 4. */
(function() {
    const LS_SEEN = 'molaop_tour_seen', LS_STEP = 'molaop_tour_step';
    function getSeen() { try { return localStorage.getItem(LS_SEEN) === 'true'; } catch(e) { return false; } }
    function setSeen(v) { try { localStorage.setItem(LS_SEEN, v); } catch(e) {} }
    function getStep() { try { return localStorage.getItem(LS_STEP); } catch(e) { return null; } }
    function setStep(v) { try { localStorage.setItem(LS_STEP, String(v)); } catch(e) {} }
    function clearStep() { try { localStorage.removeItem(LS_STEP); } catch(e) {} }

    /* TOUR_STEPS — selectors use [data-tour-target=X] per D-06/D-07.
       Copy locked verbatim from UI-SPEC §Copywriting Contract. */
    const TOUR_STEPS = [
        { id:'upload',        page:'/',        selector:'[data-tour-target="upload-area"]',
          title:'Start with your data — or a demo',
          body:"Upload your own CSV / TSV / TXT, or pick a demo dataset. For this tour we've already preloaded the PXR Agonist demo.",
          position:'right',   on_enter:null,   on_next:'startDemoPreload' },
        { id:'columns',       page:'/preview', selector:'[data-tour-target="column-preview"]',
          title:'We detect your columns automatically',
          body:"We've guessed which columns hold gene IDs, log2FC, and p-values. You can override any of them.",
          position:'top',     on_enter:null,   on_next:null },
        { id:'thresholds',    page:'/preview', selector:'[data-tour-target="thresholds"]',
          title:'Set significance thresholds',
          body:"Genes above these cutoffs count as 'significantly changed' in the over-representation test.",
          position:'right',   on_enter:null,   on_next:null },
        { id:'aop-picker',    page:'/preview', selector:'[data-tour-target="aop-picker"]',
          title:'Pick an AOP to analyse against',
          body:"We've preselected the PXR-anchored liver steatosis AOP. Click Analyse when you're ready — the tour resumes on the results page.",
          position:'bottom',  on_enter:'preselectPxrAop', on_next:null },
        { id:'results-table', page:'/analyze', selector:'[data-tour-target="results-table"]',
          title:'Enriched Key Events',
          body:'Each row is a Key Event — adjusted p-value < 0.05 means your data significantly activates it.',
          position:'top',     on_enter:null,   on_next:null },
        { id:'cytoscape',     page:'/analyze', selector:'[data-tour-target="cytoscape"]',
          title:'Visualise the AOP',
          body:'Click a KE node to see its genes. Drag to pan, scroll to zoom. Export PNG or download the network from the buttons above.',
          position:'left',    on_enter:null,   on_next:null }
    ];

    /* preselectPxrAop — dual-write typeahead (#aop_search visible + #aop_selection hidden).
       Values from config.py CASE_STUDY_AOPS["steatosis"] (RESEARCH §Pitfall 2). */
    function preselectPxrAop() {
        const sel = document.getElementById('aop_selection');
        const search = document.getElementById('aop_search');
        if (!sel || !search) return;
        sel.value = 'AOP:DEMO';
        search.value = 'PXR activation leading to liver steatosis';
        search.dispatchEvent(new Event('input', { bubbles: true }));
        search.dispatchEvent(new Event('change', { bubbles: true }));
    }

    /* startDemoPreload — step 0 on_next: POST /preview with demo_file.
       GSE90122_TO90137.tsv = "PXR Agonist 1" per config.py:34 (RESEARCH §Assumptions A4).
       columns_confirmed=true triggers auto-detected column confirmation server-side so the
       /preview response renders the column-preview, thresholds, AND aop-picker anchors
       (the latter two only mount inside `{% if volcano_data %}`). Sets step=1 (column-preview
       index) so on arrival the resume block in DOMContentLoaded fires showStep(1). */
    function startDemoPreload() {
        const form = document.createElement('form');
        form.method = 'POST'; form.action = '/preview'; form.style.display = 'none';
        [['demo_file','GSE90122_TO90137.tsv'],
         ['recommended_aops','AOP:DEMO'],
         ['columns_confirmed','true'],
         ['csrf_token',(document.querySelector('meta[name="csrf-token"]')||{}).content||'']
        ].forEach(function(pair) {
            const inp = document.createElement('input');
            inp.type = 'hidden'; inp.name = pair[0]; inp.value = pair[1];
            form.appendChild(inp);
        });
        setStep('1');
        document.body.appendChild(form);
        form.submit();
    }

    let stepIdx = -1, tipEl = null, hlEl = null, trapH = null, escH = null;

    function removeTip() {
        if (tipEl && tipEl.parentNode) tipEl.parentNode.removeChild(tipEl);
        tipEl = null;
        if (hlEl) { hlEl.classList.remove('tour-highlight'); hlEl.removeAttribute('aria-describedby'); hlEl = null; }
        if (trapH) { document.removeEventListener('keydown', trapH); trapH = null; }
        if (escH)  { document.removeEventListener('keydown', escH);  escH  = null; }
    }

    function posTip(targetEl, position, tip) {
        const r = targetEl.getBoundingClientRect();
        const sx = window.scrollX || 0, sy = window.scrollY || 0;
        const tw = 372, th = 180, g = 12;
        let top, left;
        if (position === 'right')  { top = r.top+sy+r.height/2-th/2; left = r.right+sx+g; }
        else if (position === 'left')   { top = r.top+sy+r.height/2-th/2; left = r.left+sx-tw-g; }
        else if (position === 'bottom') { top = r.bottom+sy+g;             left = r.left+sx+r.width/2-tw/2; }
        else                            { top = r.top+sy-th-g;             left = r.left+sx+r.width/2-tw/2; }
        const vpW = document.documentElement.clientWidth;
        if (left+tw > vpW+sx-8) left = vpW+sx-tw-8;
        if (left < sx+8) left = sx+8;
        if (top < sy+8) top = sy+8;
        tip.style.top = top+'px'; tip.style.left = left+'px';
        if (hlEl && hlEl !== targetEl) { hlEl.classList.remove('tour-highlight'); hlEl.removeAttribute('aria-describedby'); }
        targetEl.classList.add('tour-highlight');
        targetEl.setAttribute('aria-describedby', 'tour-tooltip-body');
        hlEl = targetEl;
    }

    function esc(s) {
        return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
    }

    function showStep(index) {
        removeTip();
        if (index < 0 || index >= TOUR_STEPS.length) return;
        const step = TOUR_STEPS[index];
        stepIdx = index; setStep(String(index));
        const targetEl = document.querySelector(step.selector);
        if (!targetEl) return;
        if (step.on_enter === 'preselectPxrAop') preselectPxrAop();
        const isLast = index === TOUR_STEPS.length - 1;
        const tip = document.createElement('div');
        tip.className = 'tour-tooltip';
        tip.setAttribute('role','dialog');
        tip.setAttribute('aria-labelledby','tour-tooltip-title');
        tip.setAttribute('aria-describedby','tour-tooltip-body');
        tip.setAttribute('data-position', step.position);
        let dots = '';
        for (let i = 0; i < TOUR_STEPS.length; i++) {
            dots += '<span class="tour-progress-dot'+(i===index?' tour-progress-dot--active':'')+'"'+
                    (i===index?' aria-current="step"':'')+' aria-hidden="true"></span>';
        }
        tip.innerHTML = '<h3 id="tour-tooltip-title" class="tour-tooltip__title">'+esc(step.title)+'</h3>'+
            '<p id="tour-tooltip-body" class="tour-tooltip__body">'+esc(step.body)+'</p>'+
            '<div class="tour-tooltip__footer"><div class="tour-progress">'+dots+'</div>'+
            '<div class="tour-actions">'+
            '<button type="button" class="btn btn--secondary" data-tour-action="skip">Skip</button>'+
            '<button type="button" class="btn btn--accent" data-tour-action="next">'+(isLast?'Finish':'Next')+'</button>'+
            '</div></div>';
        document.body.appendChild(tip); tipEl = tip;
        posTip(targetEl, step.position, tip);
        const skipBtn = tip.querySelector('[data-tour-action="skip"]');
        const nextBtn = tip.querySelector('[data-tour-action="next"]');
        if (skipBtn) skipBtn.focus();
        trapH = function(e) {
            if (e.key !== 'Tab') return;
            const fc = [skipBtn, nextBtn].filter(Boolean);
            if (!fc.length) return;
            e.preventDefault();
            const ci = fc.indexOf(document.activeElement);
            fc[(ci + (e.shiftKey ? -1 : 1) + fc.length) % fc.length].focus();
        };
        escH = function(e) { if (e.key === 'Escape') showSkipPrompt(); };
        document.addEventListener('keydown', trapH);
        document.addEventListener('keydown', escH);
        skipBtn.addEventListener('click', showSkipPrompt);
        nextBtn.addEventListener('click', handleNext);
    }

    function showSkipPrompt() {
        if (!tipEl) return;
        tipEl.innerHTML = '<p class="tour-tooltip__body">Skip this tour?</p>'+
            '<div class="tour-tooltip__footer"><div class="tour-actions">'+
            '<button type="button" class="btn btn--secondary" data-tour-action="skip-once">Skip this time</button>'+
            '<button type="button" class="btn btn--accent" data-tour-action="dont-show">Don\'t show again</button>'+
            '</div></div>';
        const b1 = tipEl.querySelector('[data-tour-action="skip-once"]');
        const b2 = tipEl.querySelector('[data-tour-action="dont-show"]');
        if (b1) { b1.focus(); b1.addEventListener('click', function() { clearStep(); removeTip(); stepIdx=-1; }); }
        if (b2) { b2.addEventListener('click', function() { setSeen('true'); clearStep(); removeTip(); stepIdx=-1; hideBanner(); }); }
    }

    function handleNext() {
        const step = TOUR_STEPS[stepIdx];
        if (stepIdx === TOUR_STEPS.length - 1) { setSeen('true'); clearStep(); removeTip(); stepIdx=-1; return; }
        if (step.on_next === 'startDemoPreload') { removeTip(); startDemoPreload(); return; }
        // Step 3 (aop-picker): user clicks Analyse themselves — dismiss tooltip, keep tour resumable
        if (stepIdx === 3) { clearStep(); removeTip(); stepIdx=-1; return; }
        showStep(stepIdx + 1);
    }

    function hideBanner() {
        const b = document.querySelector('.card--banner--onboarding');
        if (b) b.style.display = 'none';
    }

    function startTour() {
        setStep('0');
        if (window.location.pathname === '/') { showStep(0); }
        else { window.location.href = '/'; }
    }

    document.addEventListener('DOMContentLoaded', function() {
        if (getSeen()) hideBanner();

        document.addEventListener('click', function(e) {
            const btn = e.target.closest('[data-action="start-tour"]');
            if (!btn) return; e.preventDefault(); startTour();
        });
        document.addEventListener('click', function(e) {
            const lnk = e.target.closest('[data-action="dismiss-tour-banner"]');
            if (!lnk) return; e.preventDefault(); setSeen('true'); hideBanner();
        });

        // Inject hidden tour field into #enrichmentForm when tour is active (D-03)
        const form = document.getElementById('enrichmentForm');
        if (form && getStep() !== null) {
            const inp = document.createElement('input');
            inp.type='hidden'; inp.name='tour'; inp.value='1';
            form.appendChild(inp);
        }

        // Results-page resume: <meta name="tour-active"> AND step >= 4 (RESEARCH §Pitfall 1)
        const meta = document.querySelector('meta[name="tour-active"]');
        const savedNum = parseInt(getStep(), 10);
        if (meta && !isNaN(savedNum) && savedNum >= 4) { showStep(4); return; }

        // Resume on / at step 0
        const rawStep = getStep();
        if (rawStep !== null && window.location.pathname === '/') {
            const n = parseInt(rawStep, 10);
            if (!isNaN(n) && n === 0) showStep(0);
        }
        // Resume on /preview at steps 1-3
        if (rawStep !== null && window.location.pathname === '/preview') {
            const n = parseInt(rawStep, 10);
            if (!isNaN(n) && n >= 1 && n <= 3) showStep(n);
        }
    });

    window.startTour = startTour;
    window.molAopTour = { start: startTour, getStep: getStep, getSeen: getSeen };
})();
