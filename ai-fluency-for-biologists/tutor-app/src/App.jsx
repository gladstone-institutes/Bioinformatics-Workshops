import React, { useState, useEffect, useRef, useCallback } from "react";
import { Send, BookOpen, GraduationCap, RefreshCw, ChevronDown, ChevronRight, ExternalLink, Loader2, X, Sprout, Lightbulb, MessageSquare, Microscope, Rocket } from "lucide-react";
const SLIDE_URLS = {
  part1: "https://gladstone-institutes.github.io/Bioinformatics-Workshops/AI_Fluency/AI_Fluency_Part_1.html",
  part2: "https://gladstone-institutes.github.io/Bioinformatics-Workshops/AI_Fluency/AI_Fluency_Part_2.html"
};
const MODULE_COLORS = {
  introduction: { badge: "bg-slate-50 text-slate-700", dot: "bg-slate-500" },
  delegation: { badge: "bg-blue-50 text-blue-700", dot: "bg-blue-500" },
  description: { badge: "bg-emerald-50 text-emerald-700", dot: "bg-emerald-500" },
  discernment: { badge: "bg-amber-50 text-amber-700", dot: "bg-amber-500" },
  diligence: { badge: "bg-rose-50 text-rose-700", dot: "bg-rose-500" },
  "wrap-up": { badge: "bg-violet-50 text-violet-700", dot: "bg-violet-500" }
};
const MODULE_LABELS = {
  introduction: "Introduction",
  delegation: "Delegation",
  description: "Description",
  discernment: "Discernment",
  diligence: "Diligence",
  "wrap-up": "Wrap-Up"
};
const SLIDE_CONCEPTS = [
  {
    id: "intro-ai-fluency", name: "What is AI Fluency?",
    description: "The ability to work with AI effectively, efficiently, ethically, and safely.",
    content: `AI fluency is the ability to work with AI systems in ways that are effective, efficient, ethical, and safe. It's not about becoming an AI expert — it's about developing practical skills to use AI as a tool in your research. Like learning to use a new instrument, you need to know its strengths and limitations.`,
    keywords: ["ai fluency", "what is", "definition", "overview", "introduction", "basics", "getting started", "practical skills"],
    slidePath: "part1", slideFragment: "#/what-is-ai-fluency", difficultyLevel: 1, module: "introduction"
  },
  {
    id: "intro-4ds", name: "The 4 Ds Framework",
    description: "Delegation, Description, Discernment, and Diligence — a practical approach to AI.",
    content: `The 4 Ds Framework — a practical approach to working with AI effectively:
1. Delegation — Deciding *what* to ask AI to do
2. Description — Crafting *how* you ask (prompting)
3. Discernment — Evaluating *whether* the output is correct
4. Diligence — Using AI *responsibly* and ethically
These four dimensions work together, not in isolation. Good delegation informs good description; good discernment requires domain expertise; diligence applies to every step.`,
    keywords: ["4 ds", "four ds", "framework", "delegation", "description", "discernment", "diligence", "overview", "approach", "structure"],
    slidePath: "part1", slideFragment: "#/the-4-ds-framework", difficultyLevel: 1, module: "introduction"
  },
  {
    id: "intro-augmentation", name: "AI as Augmentation, Not Replacement",
    description: "AI is a thinking partner — your domain knowledge is essential.",
    content: `AI tools are thinking partners, not replacements for your expertise.
- Your domain knowledge is essential — AI doesn't know your system, your data, or your scientific question
- AI is useful for tasks like brainstorming, summarizing, drafting, and troubleshooting
- You remain responsible for the scientific validity of your work
Think of AI as a very fast, very well-read research assistant who sometimes makes things up.`,
    keywords: ["augmentation", "replacement", "thinking partner", "domain knowledge", "expertise", "role of ai", "assistant", "brainstorm", "summarize", "draft", "troubleshoot"],
    slidePath: "part1", slideFragment: "#/ai-as-augmentation-not-replacement", difficultyLevel: 1, module: "introduction"
  },
  {
    id: "intro-tools", name: "Brief Overview of Tools",
    description: "Claude, ChatGPT, Gemini — all LLMs with the same core limitations.",
    content: `Popular tools include Claude (Anthropic), ChatGPT (OpenAI), and Gemini (Google).
- All are large language models (LLMs) — trained on vast text data to understand and generate language
- They share the same core limitations: hallucinations, knowledge cutoffs, and training biases
- Relative performance varies by task and changes frequently
- They do not search the internet by default (some have optional web search features)`,
    keywords: ["tools", "claude", "chatgpt", "gemini", "llm", "large language model", "compare", "which tool", "openai", "anthropic", "google", "capabilities"],
    slidePath: "part1", slideFragment: "#/brief-overview-of-tools", difficultyLevel: 1, module: "introduction"
  },
  {
    id: "del-goal-awareness", name: "Goal and Task Awareness",
    description: "Define your task clearly before opening any AI tool.",
    content: `Before opening any AI tool, clearly define what you need:
- What is the task? e.g., summarize a paper, draft a methods section, troubleshoot a protocol
- What does success look like? Define the output you want before you start
- Do you have the expertise to evaluate the result? If you can't judge whether the AI is right, you probably shouldn't delegate it
Together with Platform Awareness, this helps you decide if AI is the right tool for the job.`,
    keywords: ["goal", "task", "awareness", "define", "what to ask", "delegate", "suitable", "appropriate", "when to use", "planning", "preparation", "evaluate"],
    slidePath: "part1", slideFragment: "#/goal-and-task-awareness", difficultyLevel: 2, module: "delegation"
  },
  {
    id: "del-platform-awareness", name: "Platform Awareness",
    description: "Understanding what AI tools can and cannot do.",
    content: `Understanding what AI tools can and cannot do:
- Context window — how much text AI can process at once. Exceeding it causes the AI to "forget" earlier content
- Knowledge cutoff — AI has no knowledge after its training date. It won't know about recent papers
- Hallucination — AI can confidently state plausible but wrong things (fake citations, invented gene names)
- Context pollution — long conversations accumulate errors; start a new chat for new tasks
- Extended thinking / research modes — newer features for complex reasoning or searching external sources`,
    keywords: ["platform", "limitation", "context window", "knowledge cutoff", "hallucination", "hallucinate", "context pollution", "can't do", "limitations", "extended thinking", "research mode", "forget", "wrong", "mistake", "token"],
    slidePath: "part1", slideFragment: "#/platform-awareness", difficultyLevel: 2, module: "delegation"
  },
  {
    id: "del-choosing-tool", name: "Choosing an AI Tool",
    description: "Consider privacy policies before choosing which AI tool to use.",
    content: `Before choosing a tool, consider the privacy policy:
- What data is stored? Is it used to train future models?
- Does your institution have approved tools or data use agreements?
- Never share sensitive, unpublished, or patient-derived data without checking your institution's policies
This is especially important for biologists who work with human subjects data, proprietary sequences, or pre-publication results.`,
    keywords: ["choosing", "privacy", "policy", "data", "sensitive", "patient", "institution", "security", "confidential", "unpublished", "proprietary", "hipaa", "phi", "nih"],
    slidePath: "part1", slideFragment: "#/choosing-an-ai-tool", difficultyLevel: 2, module: "delegation"
  },
  {
    id: "del-modalities", name: "Modalities of Task Delegation",
    description: "Automation, Augmentation, and Agency — three ways to work with AI.",
    content: `Three modes of delegating tasks to AI:
- Automation — You define exactly what needs to be done, and the AI executes it.
  Example: "Rewrite this protocol section to follow the journal's formatting guidelines, do not change the content — just the formatting."
- Augmentation — You and AI collaborate iteratively as thinking partners.
  Example: Drafting and refining an abstract together over several rounds of feedback.
- Agency — You configure the AI's role and knowledge, then let it act independently.
  Example: Giving AI your experimental goals and reagent constraints, and letting it independently design and compare multiple protocol options.
In practice, Augmentation is the most productive mode for research tasks.`,
    keywords: ["modality", "modalities", "automation", "augmentation", "agency", "mode", "collaborate", "iterate", "independent", "delegate", "how to use", "approach"],
    slidePath: "part1", slideFragment: "#/modalities-of-task-delegation", difficultyLevel: 2, module: "delegation"
  },
  {
    id: "desc-overview", name: "6 Prompting Techniques",
    description: "Overview of all six techniques for writing effective prompts.",
    content: `Six techniques for effective prompts:
1. Provide Context — Give background information the AI needs
2. Show Examples — Demonstrate the format or style you want (few-shot prompting)
3. Specify Output Constraints — Define length, format, audience, and structure
4. Break Complex Tasks into Steps — Provide sequential instructions (Process Description)
5. Ask It to Think First — Use extended thinking modes for complex reasoning
6. Define the AI's Role — Set perspective (Performance Description)
These techniques map to three types of Description: Product Description (what output you want — techniques 1-3), Process Description (how AI should approach the task — technique 4), and Performance Description (how AI should behave — techniques 5-6). Combine multiple techniques for best results. A good starting combination: context + constraints + role.`,
    keywords: ["prompt", "technique", "prompting", "overview", "six", "6", "how to prompt", "write prompt", "effective prompt", "tips", "strategy", "product description", "process description", "performance description"],
    slidePath: "part1", slideFragment: "#/prompting-techniques-overview", difficultyLevel: 2, module: "description"
  },
  {
    id: "desc-context", name: "Provide Context",
    description: "Give the AI background information — field, purpose, what's already done.",
    content: `Give the AI the background information it needs to do a good job:
- What field or subfield are you working in?
- What is the purpose of this task?
- What has already been done?
BAD PROMPT: "Summarize this paper"
GOOD PROMPT: "I'm a graduate student in cardiac biology preparing for a journal club. Summarize this paper's main findings, focusing on the signaling pathways involved in cardiomyocyte differentiation."
The good prompt gives the AI your role, your field, and what to focus on — all of which shape a much more useful response.`,
    keywords: ["context", "background", "field", "purpose", "provide context", "technique 1", "information", "specific", "vague"],
    slidePath: "part1", slideFragment: "#/provide-context", difficultyLevel: 2, module: "description"
  },
  {
    id: "desc-examples", name: "Show Examples (Few-Shot)",
    description: "Demonstrate the format/style you want by providing examples.",
    content: `Show the AI what you want by providing one or more examples of the desired output. Also called "few-shot prompting" — the AI learns from the pattern you demonstrate.
Useful for defining format, style, or level of detail.
EXAMPLE PROMPT: "Summarize each paper in this format:
[Author et al., Year] — [One-sentence finding]. Method: [Key technique]. Relevance: [Why it matters to my project]."
This is especially powerful for consistent formatting across multiple items, matching a specific writing style, or when words alone can't describe the output you want.`,
    keywords: ["examples", "few-shot", "show", "format", "style", "demonstrate", "pattern", "template", "technique 2", "few shot", "one-shot"],
    slidePath: "part1", slideFragment: "#/show-examples", difficultyLevel: 3, module: "description"
  },
  {
    id: "desc-constraints", name: "Specify Output Constraints",
    description: "Define length, format, audience, and structure for the output.",
    content: `Tell the AI exactly what the output should look like:
- Length — "in 200 words," "in 3 bullet points"
- Format — "as a table," "as a numbered list," "in paragraph form"
- Audience — "for a PI," "for a first-year graduate student," "for a grant review panel"
- Structure — "organize by: (1) hypothesis, (2) methods, (3) key results"
The more specific your constraints, the closer the output will match what you need. Vague prompts get vague answers.`,
    keywords: ["constraints", "output", "length", "format", "audience", "structure", "specify", "word count", "bullet points", "table", "technique 3", "specific"],
    slidePath: "part1", slideFragment: "#/specify-output-constraints", difficultyLevel: 2, module: "description"
  },
  {
    id: "desc-good-vs-bad", name: "Bad Prompt vs. Good Prompt",
    description: "Side-by-side comparison showing the impact of prompting techniques.",
    content: `BAD PROMPT: "Summarize this paper"
GOOD PROMPT: "Summarize this cell biology paper's experimental methodology in exactly 200 words, focusing only on protein-protein interaction techniques (e.g., co-IP, pull-downs, FRET, crosslinking). Organize your summary by:
1. Primary techniques used
2. Key reagents
3. Validation methods
Write for an audience of post-docs and graduate students."
The good prompt uses context (cell biology, protein-protein interactions), output constraints (200 words, specific techniques), and structure (three-section organization) to guide the AI. This illustrates how combining techniques from the Description module produces dramatically better results.`,
    keywords: ["bad prompt", "good prompt", "example", "compare", "comparison", "before after", "improve", "better prompt", "worse prompt", "summarize paper"],
    slidePath: "part1", slideFragment: "#/bad-prompt-vs.-good-prompt", difficultyLevel: 2, module: "description"
  },
  {
    id: "desc-steps", name: "Break Complex Tasks into Steps",
    description: "Give step-by-step instructions for multi-part tasks.",
    content: `For complex tasks, give the AI a step-by-step process to follow. This is called Process Description — defining how the AI should approach the task.
EXAMPLE PROMPT: "Analyze this paper in the following steps:
Step 1: Identify the main hypothesis
Step 2: List the experimental methods used
Step 3: Summarize the key findings
Step 4: Evaluate whether the conclusions are supported by the data"
Breaking tasks into steps prevents the AI from skipping important aspects and gives you checkpoints to evaluate each part independently.`,
    keywords: ["steps", "step by step", "break down", "complex", "sequential", "process description", "technique 4", "multi-step", "instructions", "workflow", "analyze"],
    slidePath: "part1", slideFragment: "#/break-complex-tasks-into-steps", difficultyLevel: 3, module: "description"
  },
  {
    id: "desc-think", name: "Let It Think Deeply",
    description: "Use extended thinking/reasoning mode for complex questions.",
    content: `Most AI tools now have a thinking or reasoning mode — use it for complex or multi-part questions.
- Don't need to say "think step by step" — modern tools have this built in
- Do toggle on extended thinking for harder questions; skip it for simple tasks
EXAMPLE: Turn on extended thinking mode, then ask:
"Here is a paper on CRISPR-based gene therapy for sickle cell disease. What are the key methodological limitations, and how might they affect the authors' conclusions?"
Extended thinking is most useful for: evaluating experimental design, comparing complex alternatives, identifying logical gaps, and multi-factor analysis.`,
    keywords: ["think", "thinking", "reasoning", "extended thinking", "deep", "complex", "multi-step", "technique 5", "reasoning mode", "step by step"],
    slidePath: "part1", slideFragment: "#/let-it-think-deeply", difficultyLevel: 3, module: "description"
  },
  {
    id: "desc-role", name: "Define the AI's Role",
    description: "Set the AI's perspective to shape its responses.",
    content: `Set the AI's perspective to shape its responses. This is Performance Description — defining how the AI should behave during your interaction.
EXAMPLES:
- "Explain this as a journal club presenter would to a mixed audience of biologists"
- "Write as if you are teaching a graduate student new to this subfield"
- "Act as a critical reviewer evaluating this manuscript for a high-impact journal"
Different roles produce very different outputs from the same information. A "supportive teacher" gives encouraging explanations; a "critical reviewer" points out weaknesses. Choose the role that matches what you actually need.`,
    keywords: ["role", "perspective", "act as", "persona", "performance description", "technique 6", "reviewer", "teacher", "presenter", "journal club", "tone", "style"],
    slidePath: "part1", slideFragment: "#/define-the-ais-role", difficultyLevel: 3, module: "description"
  },
  {
    id: "desc-meta", name: "Meta-Prompting",
    description: "Ask AI to help you write a better prompt.",
    content: `You can ask the AI to help you write a better prompt — this is called meta-prompting.
EXAMPLE PROMPT: "I'm trying to analyze this paper's methodology for my literature review. Can you help me craft a better prompt to get a useful summary?"
Meta-prompting is especially useful when:
- You're not sure how to structure your request
- You want to explore what's possible before committing to a specific approach
- Your initial prompt isn't giving you good results
This is a bonus technique — it's not one of the core 6, but it's a powerful strategy for improving your prompting over time.`,
    keywords: ["meta-prompting", "meta", "help write prompt", "improve prompt", "prompt engineering", "stuck", "not working", "better results", "refine", "iterate", "bonus"],
    slidePath: "part1", slideFragment: "#/bonus-ask-ai-to-help-with-your-prompt", difficultyLevel: 3, module: "description"
  },
  {
    id: "disc-types", name: "Three Types of Discernment",
    description: "Product, Process, and Performance — three lenses for evaluating AI.",
    content: `Three types of discernment for evaluating AI:
1. Product Discernment — Evaluating the quality of what AI produces
2. Process Discernment — Evaluating how the AI arrived at its output (the conversation quality)
3. Performance Discernment — Evaluating how well the AI is working for you (systematic verification)
Think of it as: What did it produce? How did it get there? How well did it work with you?`,
    keywords: ["discernment", "types", "product", "process", "performance", "evaluate", "evaluation", "assess", "judge", "quality", "overview"],
    slidePath: "part2", slideFragment: "#/three-types-of-discernment", difficultyLevel: 2, module: "discernment"
  },
  {
    id: "disc-product", name: "Product Discernment",
    description: "Evaluate AI output for accuracy, precision, completeness, and appropriateness.",
    content: `Evaluate every AI output against four criteria:
- Accuracy — Are the facts correct? Are citations real? Are gene names and pathways right?
- Precision — Is the language specific enough, or vague and generic?
- Completeness — Are key findings, methods, or caveats missing?
- Appropriateness — Is the tone, depth, and format right for your audience and purpose?
Your domain expertise is what makes this evaluation possible. A non-expert might not catch that a citation is fabricated or that a mechanism has been oversimplified. This is why you should only delegate tasks you can evaluate.`,
    keywords: ["product discernment", "accuracy", "precision", "completeness", "appropriateness", "evaluate output", "check", "verify", "quality", "correct", "wrong", "criteria"],
    slidePath: "part2", slideFragment: "#/product-discernment", difficultyLevel: 2, module: "discernment"
  },
  {
    id: "disc-errors", name: "Common AI Errors in Science",
    description: "Hallucinated citations, oversimplified mechanisms, fabricated details, and more.",
    content: `Watch out for these frequent AI mistakes in scientific contexts:
- Hallucinated citations — AI invents plausible-sounding references that don't exist (e.g., "Smith et al., 2023, Nature" when no such paper exists)
- Oversimplified mechanisms — Complex signaling cascades reduced to linear pathways
- Confused methodology — Mixing up techniques (e.g., describing Western blot results as flow cytometry data)
- Overstated conclusions — Presenting preliminary findings as established facts
- Fabricated details — Inventing sample sizes, p-values, or experimental conditions that weren't in the source
These errors are especially dangerous because they look plausible to non-experts. Your domain knowledge is the primary defense.`,
    keywords: ["errors", "mistakes", "hallucinate", "hallucination", "citation", "fabricated", "wrong", "oversimplified", "confused", "overstated", "fake", "invented", "common errors", "science", "scientific"],
    slidePath: "part2", slideFragment: "#/common-ai-errors-in-science", difficultyLevel: 3, module: "discernment"
  },
  {
    id: "disc-process", name: "Process Discernment",
    description: "Evaluate how well the AI conversation is going — not just the output.",
    content: `Evaluate how the AI is working with you, not just what it produces.
BAD SIGNS (consider stopping or starting over):
- AI consistently misunderstands technical terms in your field
- Outputs require more correction than writing from scratch would
- AI introduces new errors while fixing old ones (error propagation)
GOOD SIGNS (the collaboration is productive):
- AI helps structure your thinking in useful ways
- Iterations genuinely improve the output
- AI identifies aspects you hadn't considered
If the process isn't working, the product won't be good either. Sometimes starting a new chat with a better prompt is more efficient than continuing to fix a broken conversation.`,
    keywords: ["process discernment", "collaboration", "conversation", "working together", "productive", "iterating", "improving", "errors", "misunderstands", "starting over", "new chat"],
    slidePath: "part2", slideFragment: "#/process-discernment", difficultyLevel: 3, module: "discernment"
  },
  {
    id: "disc-when-stop", name: "When to Stop Using AI",
    description: "Recognize when AI is costing you more time than it saves.",
    content: `Consider stopping when:
- You've spent more time fixing AI output than it would take to do the task yourself
- The AI keeps making the same type of error despite reprompting
- The task requires judgment or nuance that the AI can't provide
- You find yourself accepting AI output without critically evaluating it (this is the most dangerous sign)
The goal is to save time and improve quality. If AI isn't doing both, switch approaches. This might mean trying a different prompt, a different tool, or doing the task manually.`,
    keywords: ["stop", "when to stop", "not working", "wasting time", "diminishing returns", "give up", "switch", "manual", "not helpful", "accepting without checking"],
    slidePath: "part2", slideFragment: "#/when-to-stop-using-ai-for-a-task", difficultyLevel: 3, module: "discernment"
  },
  {
    id: "disc-verification", name: "Verification Workflow",
    description: "A 5-step systematic process for checking every AI output.",
    content: `A systematic process for checking AI output — make this a habit for every AI-generated output you plan to use:
1. Compare to source — Does the AI output accurately reflect the original text?
2. Verify factual claims — Are specific statements, numbers, and names correct?
3. Check citations — Do referenced papers actually exist? Do they say what the AI claims?
4. Validate data — Are numerical values, statistics, and units accurate?
5. Assess overall quality — Is the output useful, or does it need significant revision?
This isn't about being paranoid — it's about building a reliable workflow. Just as you wouldn't submit an experiment without controls, don't submit AI-assisted work without verification.`,
    keywords: ["verification", "verify", "check", "workflow", "systematic", "citations", "facts", "data", "quality", "habit", "process", "5 steps", "five steps"],
    slidePath: "part2", slideFragment: "#/performance-discernment-verification", difficultyLevel: 3, module: "discernment"
  },
  {
    id: "dil-creation", name: "Creation Diligence",
    description: "Consider copyright, plagiarism, integrity, and data privacy.",
    content: `Be thoughtful about which AI systems you use and how you interact with them:
- Copyright and IP — AI-generated text may inadvertently reproduce copyrighted material
- Plagiarism — Submitting AI-generated text as your own without disclosure can violate academic integrity policies
- Scientific integrity — AI should support your analysis, not replace your scientific judgment
- Data privacy — Never paste unpublished data, patient information, or proprietary sequences into AI tools without checking your institution's data use policies
This is about responsible use at the point of creation — before you even start evaluating the output.`,
    keywords: ["creation diligence", "copyright", "plagiarism", "integrity", "privacy", "data", "ip", "intellectual property", "academic integrity", "sensitive", "patient", "proprietary", "responsible"],
    slidePath: "part2", slideFragment: "#/creation-diligence", difficultyLevel: 2, module: "diligence"
  },
  {
    id: "dil-transparency", name: "Transparency Diligence",
    description: "Disclose AI use in manuscripts, grants, and with collaborators.",
    content: `Be honest about AI's role in your work with everyone who needs to know:
- In manuscripts — Disclose AI use in Methods, Acknowledgments, or Author Contributions as required by the journal
- In grant proposals — AI use in grants is generally discouraged; check funder guidelines (e.g., NIH, NSF)
- With collaborators and supervisors — Let your team know when and how you used AI
- General principle: If someone would want to know AI was involved, disclose it
Transparency protects you and builds trust. It's better to over-disclose than to have someone discover undisclosed AI use later.`,
    keywords: ["transparency", "disclose", "disclosure", "manuscript", "grant", "collaborator", "supervisor", "acknowledge", "methods section", "nih", "nsf", "funder", "honest"],
    slidePath: "part2", slideFragment: "#/transparency-diligence", difficultyLevel: 2, module: "diligence"
  },
  {
    id: "dil-deployment", name: "Deployment Diligence",
    description: "A 5-step workflow: Generate, Verify, Edit, Review, Document.",
    content: `Take responsibility for verifying and vouching for every output you use or share. Follow this 5-step workflow:
1. Generate — Use AI to produce a draft or analysis
2. Verify — Check output against primary sources and your domain knowledge
3. Edit and refine — Revise for accuracy, tone, and completeness (don't just accept the draft)
4. Final human review — Read it as if you wrote it yourself. Would you stand behind it? Would you stake your reputation on it?
5. Document AI's role — Note what AI contributed for transparency
Principle: You vouch for everything you submit. If you wouldn't put your name on it without AI, don't put your name on it with AI.`,
    keywords: ["deployment", "workflow", "generate", "verify", "edit", "review", "document", "submit", "responsible", "vouch", "stand behind", "reputation", "5 steps"],
    slidePath: "part2", slideFragment: "#/deployment-diligence", difficultyLevel: 3, module: "diligence"
  },
  {
    id: "dil-journals", name: "Journal AI Policies",
    description: "Most major journals now require AI use disclosure.",
    content: `Most major journals now have AI use policies. Common themes:
- Disclosure required — Authors must declare AI-assisted writing or analysis
- AI cannot be an author — AI tools do not meet authorship criteria (no accountability)
- Human responsibility — Authors are fully responsible for all content, including AI-generated portions
- No AI-generated figures — Many journals prohibit AI-generated images (risk of misrepresentation)
- Evolving landscape — Policies change rapidly; always check your target journal's current guidelines before submission
When in doubt, check the journal's author guidelines — they're usually in the "Instructions for Authors" section.`,
    keywords: ["journal", "policy", "policies", "publication", "author", "authorship", "disclosure", "figures", "submission", "guidelines", "publisher", "evolving"],
    slidePath: "part2", slideFragment: "#/journal-ai-policies", difficultyLevel: 3, module: "diligence"
  },
  {
    id: "wrap-takeaways", name: "Key Takeaways",
    description: "The 4 Ds in summary — AI cannot replace your scientific judgment.",
    content: `The 4 Ds in summary:
1. Delegate thoughtfully — Define your task, know your platform, choose the right mode
2. Describe clearly — Use prompting techniques to get the output you need
3. Discern critically — Evaluate what AI produces, how it got there, and how it works with you
4. Act with Diligence — Use AI responsibly, transparently, and always verify before you share
AI is a powerful tool for research — but it cannot replace your judgment.`,
    keywords: ["takeaways", "summary", "key points", "recap", "conclusion", "4 ds", "four ds"],
    slidePath: "part2", slideFragment: "#/key-takeaways", difficultyLevel: 1, module: "wrap-up"
  }
];
const CONCEPTS = SLIDE_CONCEPTS.map(c => ({
  ...c,
  fullUrl: `${SLIDE_URLS[c.slidePath]}${c.slideFragment}`
}));
const LEVEL_DESCRIPTIONS = {
  1: "complete beginner who has never used AI tools",
  2: "someone who has tried AI a few times but doesn't use it regularly",
  3: "a regular AI user comfortable with basic prompting",
  4: "a daily AI user who integrates it into their workflow",
  5: "an advanced user with deep prompting and evaluation skills"
};
const SELF_REPORTED_LEVELS = [
  { level: 1, label: "I've never used AI tools", sublabel: "Complete beginner", icon: Sprout, color: "text-emerald-500", bg: "bg-emerald-50" },
  { level: 2, label: "I've tried AI a few times", sublabel: "Occasional user", icon: Lightbulb, color: "text-blue-500", bg: "bg-blue-50" },
  { level: 3, label: "I use AI regularly", sublabel: "Regular user", icon: MessageSquare, color: "text-violet-500", bg: "bg-violet-50" },
  { level: 4, label: "I use AI daily in my work", sublabel: "Daily user", icon: Microscope, color: "text-amber-500", bg: "bg-amber-50" },
  { level: 5, label: "I integrate AI deeply into my work", sublabel: "Power user", icon: Rocket, color: "text-rose-500", bg: "bg-rose-50" }
];

function retrieveSlides(question) {
  const q = question.toLowerCase();
  const qWords = q.split(/\W+/).filter(w => w.length > 2);
  const scored = CONCEPTS.map(slide => {
    let score = 0;
    for (const kw of slide.keywords) {
      if (q.includes(kw)) {
        score += kw.split(/\s+/).length * 3;
      }
    }
    const kwText = slide.keywords.join(' ');
    for (const word of qWords) {
      if (kwText.includes(word)) score += 1;
    }
    const nameLC = slide.name.toLowerCase();
    for (const word of qWords) {
      if (nameLC.includes(word)) score += 2;
    }
    return { slide, score };
  });
  scored.sort((a, b) => b.score - a.score);
  const threshold = Math.max(scored[0]?.score * 0.3, 2);
  const results = scored.filter(s => s.score >= threshold).slice(0, 4);
  if (results.length < 2) {
    return scored.slice(0, 2).filter(s => s.score > 0).map(s => s.slide);
  }
  return results.map(s => s.slide);
}

const WORKSHOP_KNOWLEDGE = `
WORKSHOP: "AI Fluency for Biologists" — Gladstone Bioinformatics Core
FRAMEWORK: The "4 Ds" of AI Fluency
Each D builds on the others. Good delegation requires understanding limitations (platform awareness). Good description requires knowing what success looks like (goal awareness). Good discernment requires domain expertise. Diligence applies at every step.
MODULE 1 — DELEGATION (Part 1)
Core question: "What should I ask AI to do?"
Key concepts:
- Goal & Task Awareness: Define task, success criteria, and your ability to evaluate before opening AI
- Platform Awareness: Context windows, knowledge cutoffs, hallucination, context pollution, extended thinking modes
- Choosing a Tool: Privacy policies, institutional data use agreements, never share sensitive/patient/unpublished data
- Three Modalities: Automation (you define, AI executes), Augmentation (collaborative iteration — most productive), Agency (AI acts independently)
MODULE 2 — DESCRIPTION (Part 1)
Core question: "How do I ask AI effectively?"
Six prompting techniques:
1. Provide Context — field, purpose, prior work. BAD: "Summarize this paper." GOOD: "I'm a graduate student in cardiac biology preparing for journal club. Summarize focusing on signaling pathways in cardiomyocyte differentiation."
2. Show Examples (Few-Shot) — demonstrate format. Example template: "[Author, Year] — [Finding]. Method: [Technique]. Relevance: [Why it matters]."
3. Specify Output Constraints — length, format, audience, structure. "200 words, as a table, for post-docs, organized by hypothesis/methods/results"
4. Break into Steps — Process Description. "Step 1: Identify hypothesis. Step 2: List methods. Step 3: Summarize findings. Step 4: Evaluate if conclusions are supported."
5. Let It Think — toggle extended thinking for complex questions (don't say "think step by step")
6. Define Role — Performance Description. "Act as a critical reviewer" vs "Act as a teacher" produces very different outputs
Bonus: Meta-prompting — ask AI to help you write a better prompt
MODULE 3 — DISCERNMENT (Part 2)
Core question: "Is the AI output correct and useful?"
Three types:
- Product Discernment: Accuracy, Precision, Completeness, Appropriateness
- Process Discernment: Is the collaboration productive? Bad: misunderstands terms, introduces new errors. Good: structures thinking, identifies new aspects
- Performance Discernment: Systematic verification workflow
Common AI errors in science: hallucinated citations, oversimplified mechanisms, confused methodology (e.g., Western blot described as flow cytometry), overstated conclusions, fabricated details (fake p-values, sample sizes)
When to stop: more time fixing than doing manually, same errors recur, accepting without evaluating
Verification workflow: 1) Compare to source, 2) Verify claims, 3) Check citations, 4) Validate data, 5) Assess quality
MODULE 4 — DILIGENCE (Part 2)
Core question: "Am I using AI responsibly?"
- Creation Diligence: copyright, plagiarism, scientific integrity, data privacy
- Transparency Diligence: disclose in manuscripts (Methods/Acknowledgments), grants (check funder), with collaborators. Rule: if someone would want to know, disclose it
- Deployment Workflow: Generate → Verify → Edit → Review ("would you stake your reputation?") → Document AI's role
- Journal Policies: disclosure required, AI can't be author, human responsibility, no AI-generated figures, policies evolving
RESOURCES:
- Anthropic AI Fluency Curriculum: anthropic.com/education (CC BY-NC-SA 4.0)
- Prompt Engineering Guide: docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/overview
`;

async function callAPI(messages, system) {
  const res = await fetch("https://api.anthropic.com/v1/messages", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ model: "claude-sonnet-4-20250514", max_tokens: 1500, system, messages })
  });
  if (!res.ok) throw new Error(`API error: ${res.status}`);
  return res.json();
}

function buildSystemPrompt(level, slides) {
  const slideDetails = slides.length > 0
    ? slides.map(s =>
      `📖 **${s.name}** [${MODULE_LABELS[s.module]}]\nLink: ${s.fullUrl}\n${s.content}`
    ).join('\n\n---\n\n')
    : "";
  const slideSection = slideDetails
    ? `\n\nMOST RELEVANT SLIDES FOR THIS QUESTION:\n(Use these as your primary source material. Teach with the examples verbatim when they fit.)\n\n${slideDetails}`
    : "";
  return `You are a friendly AI Fluency Coach for the Gladstone Bioinformatics workshop "AI Fluency for Biologists."
Student level: ${level}/5 (${LEVEL_DESCRIPTIONS[level]})
COMPLETE WORKSHOP CONTENT:
${WORKSHOP_KNOWLEDGE}
${slideSection}
HOW TO USE THE SLIDE CONTENT:
- The "MOST RELEVANT SLIDES" section above contains the specific slides most related to the student's question, with full teaching content including examples from the actual workshop
- Use the verbatim examples from slides (the BAD PROMPT / GOOD PROMPT comparisons, the biology scenarios) — they are carefully crafted teaching material
- Always link to relevant slides using the format: 📖 [Slide Name](url)
- When a slide has a concrete example that answers the student's question, lead with that example
- The "COMPLETE WORKSHOP CONTENT" section gives you the full workshop context so you can make connections across modules even when specific slides aren't retrieved
TEACHING APPROACH BY LEVEL:
- Level 1-2: Use analogies ("think of AI like a very fast research assistant who sometimes makes things up"), avoid jargon, give concrete biology examples, focus on one concept at a time
- Level 3: Balance theory with practical application, introduce technique combinations, connect concepts across modules
- Level 4-5: Be concise, focus on nuance and edge cases, discuss trade-offs, acknowledge evolving best practices
RESPONSE RULES:
1. Ground your answers in the workshop content above. Don't invent new frameworks — teach the 4 Ds
2. Include 📖 slide links for the most relevant slides (1-3 per response, not every slide)
3. Use biology-relevant examples throughout: scientific writing, literature review, experimental design, protocol development, data analysis, grant writing
4. When the student asks about prompting, show them a concrete before/after example they can adapt
5. When the student asks about discernment, emphasize that domain expertise is the key — AI literacy without domain knowledge is insufficient
6. When the student asks about diligence, be practical (not preachy) and acknowledge that policies are evolving
7. Help students see connections between modules — e.g., "This relates back to Delegation because..."
8. No emojis except 📖 for slide links
9. IMPORTANT: Do NOT start responses with "Great question!", "That's a great question!", "Excellent question!", or similar generic openers. Dive directly into the answer. Vary openings naturally.
10. Keep responses focused — answer the question asked, don't try to cover everything at once`;
}

function MarkdownRenderer({ content, isUser }) {
  const renderInline = (text) => {
    const result = [];
    let i = 0;
    const regex = /(\*\*(.+?)\*\*)|(`([^`]+)`)|(\[([^\]]+)\]\(([^)]+)\))/g;
    let lastIndex = 0;
    let match;
    while ((match = regex.exec(text)) !== null) {
      if (match.index > lastIndex) {
        result.push(<span key={i++}>{text.slice(lastIndex, match.index)}</span>);
      }
      if (match[1]) {
        result.push(<strong key={i++} className="font-semibold">{match[2]}</strong>);
      } else if (match[3]) {
        result.push(<code key={i++} className={`px-1.5 py-0.5 rounded text-sm font-mono ${isUser ? 'bg-fuchsia-900/50' : 'bg-slate-100 text-fuchsia-700'}`}>{match[4]}</code>);
      } else if (match[5]) {
        result.push(<a key={i++} href={match[7]} target="_blank" rel="noopener noreferrer" className={`font-medium underline ${isUser ? 'text-pink-100 hover:text-white' : 'text-blue-600 hover:text-fuchsia-700'}`}>{match[6]}</a>);
      }
      lastIndex = regex.lastIndex;
    }
    if (lastIndex < text.length) {
      result.push(<span key={i++}>{text.slice(lastIndex)}</span>);
    }
    return result.length > 0 ? result : text;
  };
  const lines = content.split('\n');
  const elements = [];
  let i = 0;
  let inCodeBlock = false;
  let codeBlockLang = '';
  let codeBlockContent = [];
  for (let idx = 0; idx < lines.length; idx++) {
    const line = lines[idx];
    if (line.startsWith('```')) {
      if (!inCodeBlock) {
        inCodeBlock = true;
        codeBlockLang = line.slice(3).trim() || 'text';
        codeBlockContent = [];
      } else {
        elements.push(
          <div key={i++} className="my-4 rounded-xl overflow-hidden shadow-sm">
            <div className="bg-slate-800 text-slate-400 text-xs px-4 py-2 font-medium uppercase">{codeBlockLang}</div>
            <pre className="bg-slate-900 text-slate-100 p-4 overflow-x-auto text-sm"><code>{codeBlockContent.join('\n')}</code></pre>
          </div>
        );
        inCodeBlock = false;
        codeBlockLang = '';
        codeBlockContent = [];
      }
      continue;
    }
    if (inCodeBlock) {
      codeBlockContent.push(line);
      continue;
    }
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (trimmed.startsWith('#### ')) {
      elements.push(<h4 key={i++} className={`text-sm font-semibold mt-3 mb-1 ${isUser ? 'text-white' : 'text-slate-700'}`}>{renderInline(trimmed.slice(5))}</h4>);
    } else if (trimmed.startsWith('### ')) {
      elements.push(<h3 key={i++} className={`text-sm font-semibold mt-4 mb-2 ${isUser ? 'text-white' : 'text-slate-700'}`}>{renderInline(trimmed.slice(4))}</h3>);
    } else if (trimmed.startsWith('## ')) {
      elements.push(<h2 key={i++} className={`text-base font-semibold mt-5 mb-2 ${isUser ? 'text-white' : 'text-slate-800'}`}>{renderInline(trimmed.slice(3))}</h2>);
    } else if (trimmed.startsWith('# ')) {
      elements.push(<h1 key={i++} className={`text-lg font-semibold mt-4 mb-2 ${isUser ? 'text-white' : 'text-slate-800'}`}>{renderInline(trimmed.slice(2))}</h1>);
    } else if (/^\d+\.\s/.test(trimmed)) {
      const numMatch = trimmed.match(/^(\d+)\.\s+(.+)$/);
      if (numMatch) {
        elements.push(
          <div key={i++} className="flex items-start gap-2 mb-2 ml-1">
            <span className={`font-medium ${isUser ? 'text-pink-100' : 'text-fuchsia-700'}`} style={{ minWidth: '1.5rem' }}>{numMatch[1]}.</span>
            <span className="flex-1">{renderInline(numMatch[2])}</span>
          </div>
        );
      }
    } else if (/^[-*•]\s/.test(trimmed)) {
      const text = trimmed.replace(/^[-*•]\s+/, '');
      elements.push(
        <div key={i++} className="flex items-start gap-2 mb-2 ml-1">
          <span className="w-1.5 h-1.5 rounded-full bg-fuchsia-700 mt-2 flex-shrink-0" />
          <span className="flex-1">{renderInline(text)}</span>
        </div>
      );
    } else {
      elements.push(<p key={i++} className="mb-3 last:mb-0">{renderInline(trimmed)}</p>);
    }
  }
  return <div className="text-sm leading-relaxed">{elements}</div>;
}

function Card({ children, className = "", elevation = 1, onClick }) {
  const shadows = { 0: "", 1: "shadow-md", 2: "shadow-lg", 3: "shadow-xl" };
  return <div className={`bg-white rounded-2xl ${shadows[elevation]} ${onClick ? 'cursor-pointer' : ''} ${className}`} onClick={onClick}>{children}</div>;
}

function LevelSelector({ onComplete }) {
  return (
    <div className="min-h-screen bg-slate-50 flex items-center justify-center p-6">
      <Card elevation={2} className="max-w-md w-full overflow-hidden">
        <div style={{ background: 'linear-gradient(to bottom right, #9c0366, #7a024f)' }} className="px-8 py-10 text-white">
          <div className="w-16 h-16 rounded-full flex items-center justify-center mb-6" style={{ backgroundColor: 'rgba(255,255,255,0.2)' }}><GraduationCap className="w-8 h-8" /></div>
          <p className="text-pink-100 text-sm font-medium uppercase mb-1">Gladstone Bioinformatics Core</p>
          <h1 className="text-2xl font-medium text-white">AI Fluency for Biologists</h1>
          <p className="text-pink-100 mt-1 text-lg">Interactive Tutor</p>
        </div>
        <div className="p-6">
          <p className="text-slate-600 mb-6 text-center font-medium">What's your experience with AI tools?</p>
          <div className="space-y-2">
            {SELF_REPORTED_LEVELS.map(({ level, label, sublabel, icon: Icon, color, bg }) => (
              <button key={level} onClick={() => onComplete(level)} className="w-full p-4 text-left rounded-xl border-2 border-transparent hover:border-fuchsia-700/20 hover:bg-pink-50/50 transition-all flex items-center gap-4 group">
                <div className={`w-12 h-12 rounded-full ${bg} flex items-center justify-center group-hover:scale-110 transition-transform`}><Icon className={`w-6 h-6 ${color}`} /></div>
                <div className="flex-1"><div className="text-slate-800 font-medium">{label}</div><div className="text-slate-500 text-sm">{sublabel}</div></div>
                <ChevronRight className="w-5 h-5 text-slate-300 group-hover:text-fuchsia-700" />
              </button>
            ))}
          </div>
        </div>
      </Card>
    </div>
  );
}

function ChatMessage({ message, isUser }) {
  return (
    <div className={`flex ${isUser ? 'justify-end' : 'justify-start'} mb-4`}>
      <div style={{ maxWidth: '85%' }} className={`${isUser ? 'bg-purple-200 text-black rounded-3xl rounded-br-lg px-5 py-3' : 'bg-white text-slate-800 rounded-3xl rounded-bl-lg px-5 py-4 shadow-md'}`}>
        <MarkdownRenderer content={message} isUser={isUser} />
      </div>
    </div>
  );
}

function SlideCard({ slide, onClick }) {
  const colors = MODULE_COLORS[slide.module] || MODULE_COLORS.introduction;
  return (
    <Card elevation={1} onClick={onClick} className="p-4 hover:shadow-lg transition-shadow cursor-pointer group">
      <div className="flex items-start justify-between gap-3">
        <div className="flex-1 min-w-0">
          <h4 className="font-medium text-slate-800 text-sm group-hover:text-fuchsia-700 truncate">{slide.name}</h4>
          <p className="text-xs text-slate-500 mt-1.5 line-clamp-2">{slide.description}</p>
        </div>
        <div className="w-8 h-8 rounded-full bg-slate-50 group-hover:bg-pink-50 flex items-center justify-center flex-shrink-0">
          <ExternalLink className="w-4 h-4 text-slate-400 group-hover:text-fuchsia-700" />
        </div>
      </div>
      <div className="mt-3 pt-3 border-t border-slate-100 flex items-center justify-between">
        <span className={`px-3 py-1 rounded-full text-xs font-medium ${colors.badge}`}>{MODULE_LABELS[slide.module]}</span>
        <div className="flex gap-1">{[1,2,3].map(i => <div key={i} className={`w-2 h-2 rounded-full ${i <= slide.difficultyLevel ? 'bg-fuchsia-700' : 'bg-slate-200'}`} />)}</div>
      </div>
    </Card>
  );
}

function SlideSidebar({ slides, onSlideClick }) {
  const [showAll, setShowAll] = useState(false);
  const moduleOrder = ["introduction", "delegation", "description", "discernment", "diligence", "wrap-up"];
  const groupedSlides = moduleOrder.map(mod => ({
    module: mod,
    label: MODULE_LABELS[mod],
    slides: CONCEPTS.filter(c => c.module === mod)
  })).filter(g => g.slides.length > 0);
  return (
    <div className="h-full flex flex-col bg-slate-50 border-l border-slate-200">
      <div className="p-5 border-b border-slate-200 bg-white">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 rounded-full bg-pink-50 flex items-center justify-center"><BookOpen className="w-5 h-5 text-fuchsia-700" /></div>
          <div><h3 className="font-semibold text-slate-800">Relevant Slides</h3><p className="text-xs text-slate-500">From the workshop</p></div>
        </div>
      </div>
      <div className="flex-1 overflow-y-auto p-4">
        {slides.length > 0 ? (
          <div className="space-y-3">{slides.map(slide => <SlideCard key={slide.id} slide={slide} onClick={() => onSlideClick(slide.fullUrl)} />)}</div>
        ) : (
          <div className="text-center py-12 px-4">
            <div className="w-16 h-16 rounded-full bg-slate-100 flex items-center justify-center mx-auto mb-4"><BookOpen className="w-8 h-8 text-slate-300" /></div>
            <p className="text-slate-500 text-sm">Ask a question to see relevant slides</p>
          </div>
        )}
        <div className="mt-6">
          <button onClick={() => setShowAll(!showAll)} className="flex items-center gap-2 text-sm text-slate-600 hover:text-fuchsia-700 w-full px-2 py-2 rounded-lg hover:bg-white">
            {showAll ? <ChevronDown className="w-4 h-4" /> : <ChevronRight className="w-4 h-4" />}
            <span className="font-medium">Browse All Slides</span>
            <span className="text-slate-400 ml-auto">{CONCEPTS.length}</span>
          </button>
          {showAll && (
            <div className="mt-3 space-y-4">
              {groupedSlides.map(group => {
                const colors = MODULE_COLORS[group.module];
                return (
                  <div key={group.module}>
                    <h4 className="text-xs font-semibold text-slate-400 uppercase mb-2 px-2 flex items-center gap-2">
                      <span className={`w-2 h-2 rounded-full ${colors.dot}`} />
                      {group.label} ({group.slides.length})
                    </h4>
                    <div className="space-y-0.5 max-h-48 overflow-y-auto">
                      {group.slides.map(c => <button key={c.id} onClick={() => onSlideClick(c.fullUrl)} className="w-full text-left px-3 py-2 text-sm text-slate-600 hover:text-fuchsia-700 hover:bg-white rounded-lg truncate">{c.name}</button>)}
                    </div>
                  </div>
                );
              })}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

function ChatInterface({ experienceLevel, onRestart }) {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [relevantSlides, setRelevantSlides] = useState([]);
  const [history, setHistory] = useState([]);
  const messagesEndRef = useRef(null);
  const scrollToBottom = useCallback(() => { messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' }); }, []);
  useEffect(() => { scrollToBottom(); }, [messages, scrollToBottom]);
  useEffect(() => {
    setMessages([{ role: 'assistant', content: `Welcome! I'm your AI Fluency tutor for the Gladstone Bioinformatics workshop.\n\nI'll tailor explanations to level **${experienceLevel}/5** (${LEVEL_DESCRIPTIONS[experienceLevel]}).\n\nWhat would you like to learn about?\n\n- The 4 Ds framework (Delegation, Description, Discernment, Diligence)\n- Writing effective AI prompts\n- Evaluating and verifying AI outputs\n- Responsible AI use in research\n- Or any workshop topic!` }]);
  }, [experienceLevel]);
  const handleSlideClick = (url) => window.open(url, '_blank', 'noopener,noreferrer');
  const sendMessage = async () => {
    if (!input.trim() || isLoading) return;
    const userMsg = input.trim();
    setInput('');
    const newMsgs = [...messages, { role: 'user', content: userMsg }];
    setMessages(newMsgs);
    setIsLoading(true);
    try {
      const selected = retrieveSlides(userMsg);
      setRelevantSlides(selected);
      const newHistory = [...history, { role: 'user', content: userMsg }];
      const res = await callAPI(newHistory, buildSystemPrompt(experienceLevel, selected));
      const assistantMsg = res.content?.filter(b => b.type === 'text')?.map(b => b.text)?.join('\n') || "I couldn't generate a response.";
      setMessages([...newMsgs, { role: 'assistant', content: assistantMsg }]);
      setHistory([...newHistory, { role: 'assistant', content: assistantMsg }]);
    } catch (err) {
      console.error(err);
      setMessages([...newMsgs, { role: 'assistant', content: "I'm having trouble connecting. Please try again." }]);
    } finally { setIsLoading(false); }
  };
  const handleKeyDown = (e) => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); sendMessage(); } };
  const clearChat = () => { setMessages([{ role: 'assistant', content: `Chat cleared! What would you like to learn about AI fluency?` }]); setRelevantSlides([]); setHistory([]); };
  return (
    <div className="h-screen flex flex-col bg-slate-100">
      <header className="bg-white px-4 py-3 flex items-center justify-between shadow-sm z-10">
        <div className="flex items-center gap-4">
          <div className="w-11 h-11 rounded-full flex items-center justify-center shadow-md" style={{ background: 'linear-gradient(to bottom right, #9c0366, #7a024f)' }}><GraduationCap className="w-6 h-6 text-white" /></div>
          <div>
            <h1 className="font-semibold text-slate-800 text-lg">AI Fluency Tutor</h1>
            <div className="flex items-center gap-2">
              <span className="px-3 py-1 rounded-full text-xs font-medium bg-pink-50 text-fuchsia-700">Level {experienceLevel}/5</span>
              <span className="text-xs text-slate-400">AI-powered</span>
            </div>
          </div>
        </div>
        <div className="flex items-center gap-1">
          <button onClick={clearChat} title="Clear chat" className="w-10 h-10 rounded-full flex items-center justify-center text-slate-500 hover:bg-slate-100"><X className="w-5 h-5" /></button>
          <button onClick={onRestart} title="Change level" className="w-10 h-10 rounded-full flex items-center justify-center text-slate-500 hover:bg-slate-100"><RefreshCw className="w-5 h-5" /></button>
        </div>
      </header>
      <div className="flex-1 flex overflow-hidden">
        <div className="flex-1 flex flex-col min-w-0 relative">
          <div className="flex-1 overflow-y-auto p-4 pb-32">
            {messages.map((msg, idx) => <ChatMessage key={idx} message={msg.content} isUser={msg.role === 'user'} />)}
            {isLoading && (
              <div className="flex justify-start mb-4">
                <Card elevation={1} className="px-5 py-4 flex items-center gap-3">
                  <Loader2 className="w-5 h-5 text-fuchsia-700 animate-spin" />
                  <span className="text-sm text-slate-500">Thinking...</span>
                </Card>
              </div>
            )}
            <div ref={messagesEndRef} />
          </div>
          <div className="absolute bottom-0 left-0 right-0 p-4 bg-gradient-to-t from-slate-100 via-slate-100 to-transparent pt-8">
            <Card elevation={2} className="flex items-end gap-3 p-2 pl-4">
              <textarea value={input} onChange={(e) => setInput(e.target.value)} onKeyDown={handleKeyDown} placeholder="Ask about delegation, prompting, discernment, diligence, or any workshop topic..." className="flex-1 resize-none bg-transparent py-3 focus:outline-none text-slate-800 placeholder-slate-400 text-sm" rows={1} disabled={isLoading} style={{ maxHeight: '120px' }} />
              <button onClick={sendMessage} disabled={!input.trim() || isLoading} className="w-14 h-14 rounded-full bg-fuchsia-700 text-white flex items-center justify-center shadow-lg hover:shadow-xl hover:bg-fuchsia-800 disabled:opacity-50 disabled:cursor-not-allowed transition-all"><Send className="w-5 h-5" /></button>
            </Card>
          </div>
        </div>
        <div className="hidden lg:block w-80 flex-shrink-0"><SlideSidebar slides={relevantSlides} onSlideClick={handleSlideClick} /></div>
      </div>
    </div>
  );
}

export default function GladstoneAIFluencyTutor() {
  const [experienceLevel, setExperienceLevel] = useState(null);
  if (!experienceLevel) return <LevelSelector onComplete={setExperienceLevel} />;
  return <ChatInterface experienceLevel={experienceLevel} onRestart={() => setExperienceLevel(null)} />;
}