## Data Studio Overview

### What is a Component?

A component in Data Studio is a containerized script with specified input and output files along with its defined visualizations. Each component is a separate docker with only the essential library installations and a main script that can read inputs. The main script can be in R or Python. Apart from the output files, components are also responsible for writing the files required for data visualization, along with visualization parameters.

Data Studio is fully customizable as you can select the component of your choice from the predefined list, or write your own component and add it to the list. Component templates are available to help with component creation. The Python3 and R templates are available on GitHub. You can easily download the templates. Read on about the creation of the components in detail from [here](https://drive.google.com/file/d/1rYCqOKuXc9hlRj1u-tNZoYN2_IlxnjVB/view?usp=sharing).

## Data Studio Journeys

You can explore the Polly Data Studio in two ways:

**Studio Core**

Studio core allows you to build and use your own custom workflow in Data Studio. It enables you to explore and define a workflow and ultimately visualize your data on a fully customizable dashboard and report.

**Build:** It provides the flexibility to break down the goal of a workflow into the steps that should be executed by the workflow. Use any of the components from the component library for the chosen step of the analysis. Select and add different components to grow your workflow. You can then arrange the steps to complete your workflow. If you can’t find the one you are looking for, make your own custom component and host it within your workflow.

**Visualize:** Easily interact with the parameters of the selected components and customize your data visualizations through highly configurable charts like line, bar, and pie charts, area, and bubble graphs and tables, and more.

**Report:** Ultimately create a custom version of the visualization dashboard to represent the report from the workflow. Easily annotate your report, apply styles, and color themes to make your data story work.

**Studio Preset**

Once the goal and the steps of the workflow are defined in Data Studio. It can be exported as a Studio Preset that will include a predefined series of steps that are required to complete an end-to-end process.

Studio preset is present in the form of an application with the selected components that can be interacted with. It is made with the help of Studio Core for the ingestion of high throughputs of data which often requires performing the same analysis steps over and over again.

Now if you are worrying that once the workflow is exported into a preset and cannot be customized again. We got you covered there as well. In case you need to again add a library component or a custom component within the workflow, it is fairly easy to do so.

