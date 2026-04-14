window.RevealChatBubbles = function () {
  return {
    id: "RevealChatBubbles",
    init: function (deck) {

      const bubbleClasses = ['bubble-right', 'bubble-left', 'bubble-left-2', 'bubble-left-3'];

      function isBubble(el) {
        return bubbleClasses.some(cls => el.classList.contains(cls));
      }

      function isReaction(el) {
        return el.classList.contains('reaction');
      }

      function isTypingReveal(el) {
        return el.classList.contains('typing-reveal');
      }

      // Animate a typing bubble between its two states (dots ↔ text).
      // Text stays invisible (opacity:0) during the size animation and fades in
      // only after the bubble reaches full size.
      // Returns the end height (useful for scroll calculations).
      function animateBubble(bubble, toTextRevealed) {
        // Cancel any in-progress animation on this bubble
        if (bubble._animCancel) bubble._animCancel();

        const textEl = bubble.querySelector('.bubble-text');

        // Use getBoundingClientRect for sub-pixel accurate dimensions.
        // offsetWidth/offsetHeight round to integers, which causes a visible
        // snap when the inline style is cleared back to auto.
        let r = bubble.getBoundingClientRect();
        const from = { w: r.width, h: r.height };

        toTextRevealed
          ? bubble.classList.add('text-revealed')
          : bubble.classList.remove('text-revealed');

        // Keep text invisible while we measure and animate the size change
        if (toTextRevealed && textEl) {
          textEl.style.transition = 'none';
          textEl.style.opacity = '0';
        }

        r = bubble.getBoundingClientRect();
        const to = { w: r.width, h: r.height };

        if (from.w === to.w && from.h === to.h) {
          // No size change — just fade the text in immediately
          if (toTextRevealed && textEl) fadeInText(textEl);
          return to.h;
        }

        // Use max-width/max-height rather than width/height.
        // When cleared after animation, max-width reverts to the CSS value (65%) and
        // max-height reverts to none — both of which produce the same rendered size as
        // the measured `to` dimensions (since `to` was measured under those same CSS
        // constraints). This avoids the post-animation snap caused by clearing an
        // explicit width/height whose pixel value can differ from the auto-sized result.
        // Omitting overflow:hidden also keeps the bubble tail (::after) visible throughout.
        bubble.style.maxWidth = from.w + 'px';
        bubble.style.maxHeight = from.h + 'px';
        bubble.style.transition = 'none';
        void bubble.offsetWidth; // force reflow

        bubble.style.transition = 'max-width 0.25s ease, max-height 0.25s ease';
        bubble.style.maxWidth = to.w + 'px';
        bubble.style.maxHeight = to.h + 'px';

        const expectedCount = (from.w !== to.w ? 1 : 0) + (from.h !== to.h ? 1 : 0);
        let doneCount = 0;

        function onEnd(e) {
          if (e.propertyName !== 'max-width' && e.propertyName !== 'max-height') return;
          if (++doneCount < expectedCount) return;
          bubble.removeEventListener('transitionend', onEnd);
          bubble._animCancel = null;
          bubble.style.maxWidth = '';
          bubble.style.maxHeight = '';
          bubble.style.transition = '';
          if (toTextRevealed && textEl) fadeInText(textEl);
        }

        bubble.addEventListener('transitionend', onEnd);

        bubble._animCancel = () => {
          bubble.removeEventListener('transitionend', onEnd);
          bubble._animCancel = null;
          bubble.style.maxWidth = '';
          bubble.style.maxHeight = '';
          bubble.style.transition = '';
          if (textEl) {
            textEl.style.opacity = '';
            textEl.style.transition = '';
          }
        };

        return to.h;
      }

      function fadeInText(textEl) {
        void textEl.offsetWidth; // force reflow so transition fires
        textEl.style.transition = 'opacity 0.15s ease';
        textEl.style.opacity = '';  // fall back to CSS (1)
        textEl.addEventListener('transitionend', () => {
          textEl.style.transition = '';
        }, { once: true });
      }

      const buffer = 150;

      deck.on('ready', () => {
        let hasTypingBubbles = false;

        document.querySelectorAll('.chat').forEach(chat => {
          chat.addEventListener('scroll', () => {
            chat.classList.toggle('is-scrolled', chat.scrollTop > 0);
          });

          // Link each .reaction to the bubble that precedes it.
          // Reactions may be direct children (.div syntax) or wrapped in a <p> (.span syntax).
          let bubbleCounter = 0;
          const children = Array.from(chat.children);

          // Assign IDs to all direct-child bubbles first
          children.forEach(el => {
            if (isBubble(el)) el.dataset.bubbleId = `cb-${bubbleCounter++}`;
          });

          // Find all reactions anywhere inside the chat, link to preceding bubble
          chat.querySelectorAll('.reaction').forEach(reaction => {
            // Walk up to find the direct child of chat (the flow-level container)
            let flowEl = reaction;
            while (flowEl.parentElement !== chat) flowEl = flowEl.parentElement;

            // Hide the flow-level container (reaction div or its <p> wrapper)
            flowEl.style.cssText = 'height:0;overflow:hidden;margin:0!important;padding:0!important;';

            // Find the nearest preceding sibling of flowEl that is a bubble
            let preceding = flowEl.previousElementSibling;
            while (preceding && !isBubble(preceding)) preceding = preceding.previousElementSibling;
            if (!preceding) return;

            reaction.dataset.targetBubble = preceding.dataset.bubbleId;
            if (!preceding.querySelector('.bubble-reactions')) {
              const container = document.createElement('div');
              container.className = 'bubble-reactions';
              preceding.appendChild(container);
              preceding.classList.add('has-reactions');
            }
          });

          // Process typing bubbles: insert an invisible reveal-trigger fragment after each
          const slide = chat.closest('section');
          if (!slide) return;

          const typingBubbles = Array.from(chat.querySelectorAll('.typing.fragment'))
            .sort((a, b) => parseInt(a.dataset.fragmentIndex) - parseInt(b.dataset.fragmentIndex));

          typingBubbles.forEach(bubble => {
            const currentIdx = parseInt(bubble.dataset.fragmentIndex);
            if (isNaN(currentIdx)) return;
            hasTypingBubbles = true;

            // Shift all other fragments in the slide with index > currentIdx
            slide.querySelectorAll('.fragment[data-fragment-index]').forEach(frag => {
              if (frag === bubble) return;
              const idx = parseInt(frag.dataset.fragmentIndex);
              if (idx > currentIdx) frag.dataset.fragmentIndex = idx + 1;
            });

            // Wrap bubble content so we can hide/show independently
            bubble.innerHTML =
              `<div class="bubble-text">${bubble.innerHTML}</div>` +
              `<span class="typing-indicator"><span></span><span></span><span></span></span>`;
            bubble.classList.add('is-typing');

            // Insert an invisible fragment that acts as the "reveal text" trigger
            const revealFrag = document.createElement('div');
            revealFrag.className = 'fragment typing-reveal';
            revealFrag.dataset.fragmentIndex = currentIdx + 1;
            revealFrag.dataset.targetTypingBubble = bubble.dataset.bubbleId;
            revealFrag.style.cssText = 'height:0;overflow:hidden;margin:0!important;padding:0!important;';
            bubble.insertAdjacentElement('afterend', revealFrag);
          });
        });

        // Re-sync Reveal.js fragment state after DOM modifications
        if (hasTypingBubbles) deck.sync();
      });

      deck.on('fragmentshown', (event) => {
        const fragment = event.fragment;
        const chat = fragment.closest('.chat');
        if (!chat) return;

        if (isReaction(fragment)) {
          const bubble = chat.querySelector(`[data-bubble-id="${fragment.dataset.targetBubble}"]`);
          if (!bubble) return;
          const pill = document.createElement('span');
          pill.className = 'reaction-pill';
          pill.textContent = fragment.textContent.trim();
          fragment._reactionPill = pill;
          bubble.querySelector('.bubble-reactions').appendChild(pill);
          return;
        }

        if (isTypingReveal(fragment)) {
          const bubble = chat.querySelector(`[data-bubble-id="${fragment.dataset.targetTypingBubble}"]`);
          if (!bubble) return;
          const endHeight = animateBubble(bubble, true);
          // Scroll using the known end height (bubble.offsetHeight is mid-animation)
          const fragBottom = bubble.offsetTop + endHeight;
          const visibleBottom = chat.scrollTop + chat.clientHeight - buffer;
          if (fragBottom > visibleBottom) {
            chat.scrollTo({ top: fragBottom - chat.clientHeight + buffer, behavior: 'smooth' });
          }
          return;
        }

        if (!isBubble(fragment)) return;
        const fragBottom = fragment.offsetTop + fragment.offsetHeight;
        const visibleBottom = chat.scrollTop + chat.clientHeight - buffer;
        if (fragBottom > visibleBottom) {
          chat.scrollTo({ top: fragBottom - chat.clientHeight + buffer, behavior: 'smooth' });
        }
      });

      deck.on('fragmenthidden', (event) => {
        const fragment = event.fragment;
        const chat = fragment.closest('.chat');
        if (!chat) return;

        if (isReaction(fragment)) {
          if (fragment._reactionPill) {
            const pill = fragment._reactionPill;
            fragment._reactionPill = null;
            pill.classList.add('is-hiding');
            pill.addEventListener('animationend', () => pill.remove(), { once: true });
          }
          return;
        }

        if (isTypingReveal(fragment)) {
          const bubble = chat.querySelector(`[data-bubble-id="${fragment.dataset.targetTypingBubble}"]`);
          if (!bubble) return;
          animateBubble(bubble, false);
          return;
        }

        if (!isBubble(fragment)) return;
        const targetTop = Math.max(0, fragment.offsetTop - chat.clientHeight);
        if (chat.scrollTop > targetTop) {
          chat.scrollTo({ top: targetTop, behavior: 'smooth' });
        }
      });

    }
  };
};
